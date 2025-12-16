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
        SortSct(ctInfo,ctFolder); 
    end
end

modalitylist = {'RTDOSE','RTSTRUCT','RTPLAN'};

for modality = 1:length(modalitylist)
    ctCollection = ctInfo; 
    rowIndices = strcmp(ctCollection.Modality, modalitylist(modality)); 
    selectedfiles = vertcat(ctCollection.Filenames{rowIndices}); 
    trimmedCollection = ctCollection(rowIndices, :); 
    rt_file = trimmedCollection.Filenames{1}{1}; 
    finaldirectory = fullfile(ctFolder,'sct'); 
    
    copyfile(rt_file,finaldirectory);  

end
%% Functions

function void = SortSct(ctInfo,ctFolder)

    % SCT
    parent = ctFolder;
    sctDir = fullfile(parent,'sct');
    mkdir(sctDir);
    rowIndex = strcmp(ctInfo.SeriesDescription, 'sct'); 
    sctInfo = ctInfo(rowIndex,:); 
    
    sctFiles = sctInfo.Filenames{1}; 
    for k = 1:length(sctFiles);
        file = sctFiles(k);
        [~, name, ext] = fileparts(file);
        newFile = fullfile(sctDir,name +ext);
        movefile(file,newFile);
    end
end
% Other modalities
% Duplicates need to be investigated
% This can probably be simplified actually since all these files are in batches of one. Whoops

% 
% for i=1:length(modalitylist)
%     rowIndex = strcmp(ctInfo.Modality,modalitylist{i}); 
%     modalityInfo = ctInfo(rowIndex,:); 
%     for k = 1:height(modalityInfo) 
%         modalityFiles = modalityInfo.Filenames{k};
%         dirname = strcat(modalitylist{i},int2str(k));
%         modalitykdir = fullfile(parent,dirname); 
%         mkdir(modalitykdir); 
%         for j = 1:length(modalityFiles)
% 
%             file = modalityFiles(j);
%             [~, name, ext] = fileparts(file);
%             newFile = fullfile(modalitykdir,name +ext);
%             movefile(file,newFile);
%         end
%     end
% end

% Image confirmation


