classdef DoseSource
%DOSESOURCE  Where a SimResult's reconstructed dose comes from.
%
%   Drives ethos.Simulator.resolve():
%       Simulate - always run k-Wave (run_single_field_simulation), ignore any
%                  cached recon on disk.
%       Load     - require a hash-matching recon on disk (load_recon_dose_data);
%                  error if none is present.
%       Auto     - load the on-disk recon when one matching the CONFIG hash
%                  exists, otherwise fall back to simulating. Default for
%                  studies: re-runnable at zero k-Wave cost once recons exist,
%                  and self-healing when they do not.
%
%   See also: ethos.Simulator, ethos.SimConfig, ethos.StudyRunner

    enumeration
        Simulate
        Load
        Auto
    end

    methods (Static)
        function ds = fromChar(s)
        %FROMCHAR  Parse 'simulate' | 'load' | 'auto' (case-insensitive).
            if isa(s, 'ethos.DoseSource')
                ds = s;
                return;
            end
            switch lower(char(s))
                case {'simulate', 'sim'}
                    ds = ethos.DoseSource.Simulate;
                case 'load'
                    ds = ethos.DoseSource.Load;
                case 'auto'
                    ds = ethos.DoseSource.Auto;
                otherwise
                    error('ethos:DoseSource:BadValue', ...
                        'Unknown DoseSource ''%s'' (use simulate | load | auto).', ...
                        char(s));
            end
        end
    end
end
