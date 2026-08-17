function Data1024 = addMultiplexing(dataPath, numPulses)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% LOADING DATA %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(dataPath)
    % IRAI Low-frequency 2D Array
    
    for i=1:4
        PAD(:,:,i)=mean(double(PAData1(:,:,i:4:numPulses)),3);
    end


    %PAD=double(PAData(:,:,1:4));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% ADD MULTIPLEXING %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('channelswitch.mat')
    
    frame = double(PAD(1:640,:,1:4));
    
    Data1024=zeros(640,32,32);
    
    Data256=zeros(640,256,4);
    
    Data256=frame(:,po,:);
    
    Data256=reshape(Data256,[640 16 16 4]);
    
    Selection=[1 2;3 4];
    
    for sel1 = 0 : 1
        
        for sel2 = 0 : 1
            
            for i=1:16
                for j=1:16
               
                    Data1024(:,i*2 + sel1 - 1,j*2 - sel2)=Data256(:,i,j,Selection(sel1+1,sel2+1));
                end
            end
                   
        end
    end
    A=reshape(Data1024,[640 1024]);
    B=detrend(A);
    DataDetrend=reshape(B,[640 32 32] );
    disp("multiplexing finished");
end
