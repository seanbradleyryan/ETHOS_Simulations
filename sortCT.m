% Choose ethos export directory
% This could loop over directories once more data is exported

patient_ids = {'1194203'}; 
sessions = {'Session_1','Session_2'}; 
imagepath = 'Pancreas/Session_1/';

for id = 1:length(patient_ids)
    for session = 1:length(sessions)
        ctFolder = fullfile(pwd, 'EthosExports', patient_ids{id},'Pancreas',sessions{session});  
        if ~isfolder(ctFolder)
            print("Warning: Empty session")
            continue; 
        end
        ctInfo = dicomCollection(ctFolder);
        % Print diagnostic information about the collection
        fprintf('\n=== Analyzing Patient %s, %s ===\n', patient_ids{id}, sessions{session});
        PrintCollectionInfo(ctInfo);
        sctSeriesUID = SortSct(ctInfo,ctFolder); 
        SortRTFiles(ctInfo, ctFolder, sctSeriesUID);    
    end
end

function sctSeriesUID = SortSct(ctInfo, ctFolder)
    % SCT - Sort and return the SeriesInstanceUID for matching
    parent = ctFolder;
    sctDir = fullfile(parent,'sct');
    if ~isfolder(sctDir)
        mkdir(sctDir);
    end
    
    rowIndex = strcmp(ctInfo.SeriesDescription, 'sct'); 
    sctInfo = ctInfo(rowIndex,:); 
    
    if height(sctInfo) == 0
        warning('No SCT series found');
        sctSeriesUID = '';
        return;
    end
    
    % Get the SeriesInstanceUID from the first SCT file
    sctFiles = sctInfo.Filenames{1}; 
    if ~isempty(sctFiles)
        firstFile = sctFiles{1};
        sctMetadata = dicominfo(firstFile);
        sctSeriesUID = sctMetadata.SeriesInstanceUID;
        fprintf('SCT SeriesInstanceUID: %s\n', sctSeriesUID);
        fprintf('SCT Series Date/Time: %s / %s\n', ...
            sctMetadata.SeriesDate, sctMetadata.SeriesTime);
    else
        sctSeriesUID = '';
    end
    
    % Move SCT files
    for k = 1:length(sctFiles)
        file = sctFiles{k};
        [~, name, ext] = fileparts(file);
        newFile = fullfile(sctDir, strcat(name, ext));
        if exist(file, 'file')
            movefile(file, newFile);
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function void = SortRTFiles(ctInfo, ctFolder, sctSeriesUID)
    % Sort RTDOSE, RTSTRUCT, RTPLAN based on instance relationships
    
    modalitylist = {'RTSTRUCT','RTPLAN','RTDOSE'};
    
    % Extract RT file information with metadata
    rtFileInfo = struct();
    
    for modIdx = 1:length(modalitylist)
        modality = modalitylist{modIdx};
        rowIndices = strcmp(ctInfo.Modality, modality); 
        
        if ~any(rowIndices)
            fprintf('No %s files found\n', modality);
            continue;
        end
        
        selectedRows = ctInfo(rowIndices, :);
        
        % Store information for each file
        rtFileInfo.(modality) = {};
        
        for fileIdx = 1:height(selectedRows)
            fileCell = selectedRows.Filenames{fileIdx};
            if ~isempty(fileCell) && ~isempty(fileCell{1})
                filePath = fileCell{1};
                metadata = dicominfo(filePath);
                
                fileInfo = struct();
                fileInfo.FilePath = filePath;
                fileInfo.SeriesInstanceUID = metadata.SeriesInstanceUID;
                fileInfo.SeriesDate = metadata.SeriesDate;
                fileInfo.SeriesTime = metadata.SeriesTime;
                fileInfo.SeriesDescription = metadata.SeriesDescription;
                
                % Get referenced series/instance UIDs
                if strcmp(modality, 'RTSTRUCT') && isfield(metadata, 'ReferencedFrameOfReferenceSequence')
                    % RTSTRUCT references CT series
                    refSeq = metadata.ReferencedFrameOfReferenceSequence;
                    if isfield(refSeq, 'Item_1') && isfield(refSeq.Item_1, 'RTReferencedStudySequence')
                        studySeq = refSeq.Item_1.RTReferencedStudySequence;
                        if isfield(studySeq, 'Item_1') && isfield(studySeq.Item_1, 'RTReferencedSeriesSequence')
                            seriesSeq = studySeq.Item_1.RTReferencedSeriesSequence;
                            if isfield(seriesSeq, 'Item_1') && isfield(seriesSeq.Item_1, 'SeriesInstanceUID')
                                fileInfo.ReferencedSeriesUID = seriesSeq.Item_1.SeriesInstanceUID;
                            end
                        end
                    end
                elseif strcmp(modality, 'RTPLAN') && isfield(metadata, 'ReferencedStructureSetSequence')
                    % RTPLAN references RTSTRUCT
                    refSeq = metadata.ReferencedStructureSetSequence;
                    if isfield(refSeq, 'Item_1') && isfield(refSeq.Item_1, 'ReferencedSOPInstanceUID')
                        fileInfo.ReferencedSOPInstanceUID = refSeq.Item_1.ReferencedSOPInstanceUID;
                    end
                elseif strcmp(modality, 'RTDOSE') && isfield(metadata, 'ReferencedRTPlanSequence')
                    % RTDOSE references RTPLAN
                    refSeq = metadata.ReferencedRTPlanSequence;
                    if isfield(refSeq, 'Item_1') && isfield(refSeq.Item_1, 'ReferencedSOPInstanceUID')
                        fileInfo.ReferencedSOPInstanceUID = refSeq.Item_1.ReferencedSOPInstanceUID;
                    end
                end
                
                rtFileInfo.(modality){end+1} = fileInfo;
            end
        end
    end
    
    % Find the best matching set of RT files
    [selectedRTStruct, selectedRTPlan, selectedRTDose] = FindMatchingRTSet(rtFileInfo, sctSeriesUID);
    
    % Copy selected files to sct directory
    sctDir = fullfile(ctFolder, 'sct');
    if ~isfolder(sctDir)
        mkdir(sctDir);
    end
    
    if ~isempty(selectedRTStruct)
        copyfile(selectedRTStruct.FilePath, sctDir);
        fprintf('Copied RTSTRUCT: %s (Date: %s)\n', selectedRTStruct.SeriesDescription, selectedRTStruct.SeriesDate);
    end
    
    if ~isempty(selectedRTPlan)
        copyfile(selectedRTPlan.FilePath, sctDir);
        fprintf('Copied RTPLAN: %s (Date: %s)\n', selectedRTPlan.SeriesDescription, selectedRTPlan.SeriesDate);
    end
    
    if ~isempty(selectedRTDose)
        copyfile(selectedRTDose.FilePath, sctDir);
        fprintf('Copied RTDOSE: %s (Date: %s)\n', selectedRTDose.SeriesDescription, selectedRTDose.SeriesDate);
    end
end





function [selectedStruct, selectedPlan, selectedDose] = FindMatchingRTSet(rtFileInfo, sctSeriesUID)
    % Find a matching set of RTSTRUCT, RTPLAN, and RTDOSE
    % Priority: 1) Match to SCT, 2) Most recent complete set, 3) Most recent individual files
    
    selectedStruct = [];
    selectedPlan = [];
    selectedDose = [];
    
    % Get all RT structures
    rtStructs = [];
    rtPlans = [];
    rtDoses = [];
    
    if isfield(rtFileInfo, 'RTSTRUCT')
        rtStructs = rtFileInfo.RTSTRUCT;
    end
    if isfield(rtFileInfo, 'RTPLAN')
        rtPlans = rtFileInfo.RTPLAN;
    end
    if isfield(rtFileInfo, 'RTDOSE')
        rtDoses = rtFileInfo.RTDOSE;
    end
    
    % Strategy 1: Find RTSTRUCT that references the SCT
    if ~isempty(sctSeriesUID) && ~isempty(rtStructs)
        for i = 1:length(rtStructs)
            if isfield(rtStructs{i}, 'ReferencedSeriesUID') && ...
               strcmp(rtStructs{i}.ReferencedSeriesUID, sctSeriesUID)
                selectedStruct = rtStructs{i};
                fprintf('Found RTSTRUCT referencing SCT\n');
                break;
            end
        end
    end
    
    % If no match to SCT, select most recent RTSTRUCT
    if isempty(selectedStruct) && ~isempty(rtStructs)
        selectedStruct = GetMostRecent(rtStructs);
        fprintf('Selected most recent RTSTRUCT\n');
    end
    
    % Find RTPLAN that references the selected RTSTRUCT
    if ~isempty(selectedStruct) && ~isempty(rtPlans)
        structSOPInstanceUID = dicominfo(selectedStruct.FilePath).SOPInstanceUID;
        for i = 1:length(rtPlans)
            if isfield(rtPlans{i}, 'ReferencedSOPInstanceUID') && ...
               strcmp(rtPlans{i}.ReferencedSOPInstanceUID, structSOPInstanceUID)
                selectedPlan = rtPlans{i};
                fprintf('Found RTPLAN referencing selected RTSTRUCT\n');
                break;
            end
        end
    end
    
    % If no matching plan, select most recent
    if isempty(selectedPlan) && ~isempty(rtPlans)
        selectedPlan = GetMostRecent(rtPlans);
        fprintf('Selected most recent RTPLAN\n');
    end
    
    % Find RTDOSE that references the selected RTPLAN
    if ~isempty(selectedPlan) && ~isempty(rtDoses)
        planSOPInstanceUID = dicominfo(selectedPlan.FilePath).SOPInstanceUID;
        for i = 1:length(rtDoses)
            if isfield(rtDoses{i}, 'ReferencedSOPInstanceUID') && ...
               strcmp(rtDoses{i}.ReferencedSOPInstanceUID, planSOPInstanceUID)
                selectedDose = rtDoses{i};
                fprintf('Found RTDOSE referencing selected RTPLAN\n');
                break;
            end
        end
    end
    
    % If no matching dose, select most recent
    if isempty(selectedDose) && ~isempty(rtDoses)
        selectedDose = GetMostRecent(rtDoses);
        fprintf('Selected most recent RTDOSE\n');
    end
end

function mostRecent = GetMostRecent(fileList)
    % Get the most recent file based on SeriesDate and SeriesTime
    mostRecent = [];
    
    if isempty(fileList)
        return;
    end
    
    mostRecent = fileList{1};
    mostRecentDateTime = str2double([mostRecent.SeriesDate mostRecent.SeriesTime]);
    
    for i = 2:length(fileList)
        currentDateTime = str2double([fileList{i}.SeriesDate fileList{i}.SeriesTime]);
        if currentDateTime > mostRecentDateTime
            mostRecent = fileList{i};
            mostRecentDateTime = currentDateTime;
        end
    end
end

function void = PrintCollectionInfo(ctInfo)
    % Print diagnostic information about the DICOM collection
    
    fprintf('\n--- DICOM Collection Summary ---\n');
    fprintf('Total files: %d\n', height(ctInfo));
    
    % Group by modality
    modalities = unique(ctInfo.Modality);
    for i = 1:length(modalities)
        modality = modalities{i};
        count = sum(strcmp(ctInfo.Modality, modality));
        fprintf('  %s: %d series\n', modality, count);
    end
    
    % Print detailed info for RT files
    fprintf('\n--- RT File Details ---\n');
    rtModalities = {'RTSTRUCT', 'RTPLAN', 'RTDOSE'};
    
    for i = 1:length(rtModalities)
        modality = rtModalities{i};
        rowIndices = strcmp(ctInfo.Modality, modality);
        
        if any(rowIndices)
            fprintf('\n%s Files:\n', modality);
            selectedRows = ctInfo(rowIndices, :);
            
            for j = 1:height(selectedRows)
                fileCell = selectedRows.Filenames{j};
                if ~isempty(fileCell) && ~isempty(fileCell{1})
                    metadata = dicominfo(fileCell{1});
                    fprintf('  [%d] Series: %s\n', j, metadata.SeriesDescription);
                    fprintf('      Date/Time: %s / %s\n', metadata.SeriesDate, metadata.SeriesTime);
                    fprintf('      SeriesInstanceUID: %s\n', metadata.SeriesInstanceUID);
                    
                    % Print reference information if available
                    if strcmp(modality, 'RTSTRUCT') && isfield(metadata, 'ReferencedFrameOfReferenceSequence')
                        try
                            refUID = metadata.ReferencedFrameOfReferenceSequence.Item_1.RTReferencedStudySequence.Item_1.RTReferencedSeriesSequence.Item_1.SeriesInstanceUID;
                            fprintf('      References CT Series: %s\n', refUID);
                        catch
                            fprintf('      References: (structure unavailable)\n');
                        end
                    elseif strcmp(modality, 'RTPLAN') && isfield(metadata, 'ReferencedStructureSetSequence')
                        refUID = metadata.ReferencedStructureSetSequence.Item_1.ReferencedSOPInstanceUID;
                        fprintf('      References RTSTRUCT: %s\n', refUID);
                    elseif strcmp(modality, 'RTDOSE') && isfield(metadata, 'ReferencedRTPlanSequence')
                        refUID = metadata.ReferencedRTPlanSequence.Item_1.ReferencedSOPInstanceUID;
                        fprintf('      References RTPLAN: %s\n', refUID);
                    end
                end
            end
        end
    end
    
    % Print SCT info
    fprintf('\n--- SCT Series ---\n');
    sctRows = strcmp(ctInfo.SeriesDescription, 'sct');
    if any(sctRows)
        sctInfo = ctInfo(sctRows, :);
        fileCell = sctInfo.Filenames{1};
        if ~isempty(fileCell) && ~isempty(fileCell{1})
            metadata = dicominfo(fileCell{1});
            fprintf('  Series: %s\n', metadata.SeriesDescription);
            fprintf('  Date/Time: %s / %s\n', metadata.SeriesDate, metadata.SeriesTime);
            fprintf('  SeriesInstanceUID: %s\n', metadata.SeriesInstanceUID);
        end
    else
        fprintf('  No SCT series found\n');
    end
    
    fprintf('\n================================\n\n');
end
