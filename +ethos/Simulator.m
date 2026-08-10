classdef Simulator < handle
%SIMULATOR  Thin orchestrator over the ETHOS physics functions.
%
%   PURPOSE:
%   Sequences the existing pipeline functions and, crucially, owns the
%   load-or-simulate decision so the analysis layer never has to. resolve()
%   returns a SimResult whether the dose came from disk or a fresh k-Wave run,
%   driven by cfg.dose_source (Auto | Load | Simulate). No physics lives here;
%   every method delegates to run_single_field_simulation / load_recon_dose_data.
%
%   USAGE:
%       sim = ethos.Simulator(ethos.SimConfig.standalone().withOverrides( ...
%               'working_dir', WD, 'patient_id', '1194203', 'session', 'Session_1'));
%       c   = ethos.SimCase.list('1194203', 'Session_1', sim.cfg, 'Beam', 15, 'Segment', 112);
%       r   = sim.resolve(c);          % load if cached, else simulate
%       rs  = sim.resolveMany(cases);  % parfor over a case array
%
%   NOTE (noise-ensemble null): the redraw-noise ensemble in
%   study_pass_rates_individual depends on forward-bundle internals that are not
%   exposed as standalone functions, so it is not reproduced here; that study
%   still owns it. The class layer covers the deterministic forward+recon path.
%
%   See also: ethos.SimConfig, ethos.SimCase, ethos.SimResult,
%             run_single_field_simulation, load_recon_dose_data

    properties
        cfg              % ethos.SimConfig
        save_recon logical = false   % write <base>_recon_<hash>.mat after run()
    end

    methods
        function obj = Simulator(cfg)
            if nargin < 1 || isempty(cfg)
                cfg = ethos.SimConfig();
            end
            obj.cfg = cfg;
        end

        function r = resolve(obj, simCase)
        %RESOLVE  Load-or-simulate one field per cfg.dose_source. Returns a SimResult.
            ds = obj.cfg.dose_source;
            switch ds
                case ethos.DoseSource.Load
                    r = obj.loadField(simCase);
                case ethos.DoseSource.Simulate
                    r = obj.run(simCase);
                case ethos.DoseSource.Auto
                    try
                        r = obj.loadField(simCase);
                        fprintf('[Simulator] Loaded cached recon for %s.\n', simCase.label());
                    catch ME
                        fprintf(['[Simulator] No cached recon for %s (%s); ' ...
                            'simulating.\n'], simCase.label(), ME.message);
                        r = obj.run(simCase);
                    end
                otherwise
                    error('ethos:Simulator:BadSource', 'Unknown dose_source.');
            end
        end

        function results = resolveMany(obj, cases)
        %RESOLVEMANY  parfor resolve over a SimCase array -> SimResult array.
        %   Iterates value SimCases and rebuilds a Simulator per worker so no
        %   handle/mutable state crosses the parfor boundary (statelessness rule).
            % Force plotting off in workers: figures from a parfor do not render.
            cfg   = obj.cfg.withOverrides('plot_results', false); %#ok<PROPLC>
            saveR = obj.save_recon;
            n     = numel(cases);
            results = repmat(ethos.SimResult(), 1, n);
            parfor i = 1:n
                sim = ethos.Simulator(cfg);
                sim.save_recon = saveR;
                results(i) = sim.resolve(cases(i));
            end
        end

        function r = run(obj, simCase)
        %RUN  Force a k-Wave forward + reconstruction for one field.
            sc   = simCase.materialize(obj.cfg);
            cfgS = obj.cfg.toStruct();

            [recon_dose, sim_results] = run_single_field_simulation( ...
                sc.field_dose, sc.cbct_resampled, sc.medium, ...
                sc.beam_metadata, cfgS, sc.precomputed_sensor);

            r = ethos.SimResult.fromSim(sc.field_dose, sc.cbct_resampled, ...
                sc.medium, recon_dose, sim_results, obj.cfg);

            if ~isempty(sc.precomputed_sensor) && ...
                    isfield(sc.precomputed_sensor, 'sensor_mask')
                r = r.withSensor(sc.precomputed_sensor.sensor_mask, ...
                    sc.precomputed_sensor.sensor_info);
            end

            if obj.save_recon
                obj.writeRecon(sc, recon_dose);
            end
        end

        function r = runBlind(obj, fwdCase, refCase)
        %RUNBLIND  Forward on the TRUE anatomy, invert on a REFERENCE geometry.
        %   Models live IRAI: waves propagate through fwdCase's true CT, but time
        %   reversal runs on refCase's geometry (the only one we have to invert
        %   on). Wraps run_second_field_simulation (medium_recon = refCase medium).
            scF  = fwdCase.materialize(obj.cfg);
            scR  = refCase.materialize(obj.cfg);
            cfgS = obj.cfg.toStruct();

            [recon_dose, sim_results] = run_second_field_simulation( ...
                scF.field_dose, scF.cbct_resampled, scF.medium, scR.medium, ...
                scF.beam_metadata, cfgS, scF.precomputed_sensor);

            r = ethos.SimResult.fromSim(scF.field_dose, scF.cbct_resampled, ...
                scF.medium, recon_dose, sim_results, obj.cfg);
            if ~isempty(scF.precomputed_sensor) && ...
                    isfield(scF.precomputed_sensor, 'sensor_mask')
                r = r.withSensor(scF.precomputed_sensor.sensor_mask, ...
                    scF.precomputed_sensor.sensor_info);
            end
        end

        function r = loadField(obj, simCase)
        %LOADFIELD  Load one field's cached recon via load_recon_dose_data.
            if isempty(simCase.beam)
                error('ethos:Simulator:NoBeam', ...
                    'Loading requires a beam number on the SimCase.');
            end
            args = {'Mode', 'single', 'Beam', simCase.beam};
            if ~isempty(simCase.segment)
                args = [args, {'Segment', simCase.segment}];
            end
            if ~isempty(simCase.plan_type)
                args = [args, {'PlanType', simCase.plan_type}];
            end
            if ~isempty(simCase.ct_label)
                args = [args, {'CTLabel', simCase.ct_label}];
            end
            out = load_recon_dose_data(simCase.patient_id, simCase.session, ...
                obj.cfg.toStruct(), args{:});
            r = ethos.SimResult.fromField(out, 1);
        end

        function r = loadTotal(obj, patient_id, session)
        %LOADTOTAL  Load the summed total recon via load_recon_dose_data.
            out = load_recon_dose_data(patient_id, session, obj.cfg.toStruct(), ...
                'Mode', 'total');
            r = ethos.SimResult.fromTotal(out);
        end
    end

    methods (Access = private)
        function writeRecon(obj, sc, recon_dose)
        %WRITERECON  Cache the recon next to the input dose (pipeline naming).
            hash8   = obj.cfg.hash();
            sim_dir = fullfile(obj.cfg.working_dir, 'SimulationResults', ...
                sc.patient_id, sc.session, obj.cfg.gruneisen_method);
            if ~isfolder(sim_dir), mkdir(sim_dir); end
            if ~isempty(sc.source_mat_filename)
                [~, base] = fileparts(sc.source_mat_filename);
            else
                base = sprintf('dose_%s_%s_B%s_%s', sc.patient_id, sc.session, ...
                    mat2str(sc.beam), mat2str(sc.segment));
            end
            out_path    = fullfile(sim_dir, sprintf('%s_recon_%s.mat', base, hash8)); %#ok<NASGU>
            config_hash = hash8; %#ok<NASGU>
            save(out_path, 'recon_dose', 'config_hash', '-v7.3');
            fprintf('[Simulator] Saved recon: %s\n', out_path);
        end
    end
end
