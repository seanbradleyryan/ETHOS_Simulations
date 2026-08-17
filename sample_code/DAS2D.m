function DAS2D(Data1024, outfile, img_size, freq_rate)
    %  close all
    % clear all;
    MyPar = parpool; % parallel computing
    
    %load('Treatment plan 3D conformal 2.mat') %load data file
    %load("Patient2_2.mat");
    PAD = Data1024;
    
    wnd=0;
    air=dir;
    
    Dr=img_size; % image size 
    sc=1;
    fs=freq_rate;%%fs=0.7812*4 MhZ
    
    nsr=0.3;
    PAD;
    
    % C=zeros(1,14);
    % C(2)=1;
    % C(13)=-1;
    % for l=1:96
    %     PAD(:,l) = deconvwnr(real(PAD(:,l))',C,nsr)';% Wiener deconvolution 
    % end
    
    elmspace=0.85*1.54/0.35; % space between two elements
    
    
    im_pad=zeros(2*Dr+1,2*Dr+1,2*Dr+1); % image size in pixels
    
    im_pad_dump=im_pad;
    
    parfor iii=3:30
        iii
        im_pad = im_pad+pcalculation(iii,Dr,elmspace,PAD);
    
    end
    R=Dr*1.54/fs;
    xa=1:2*R/2/Dr:2*R+1;
    
    
    
    parfor i=1:(2*Dr+1)
    C1(:,i,:)=envelopeDetection(squeeze(im_pad(:,i,:)));
    end
    
    delete(MyPar)
    % scaling applied, C1 is the actual data
    figure;imagesc(2*xa,2*xa,squeeze(C1(:,:,Dr)));
    axis image
    save(outfile, 'C1', '-v7.3');
    %save('results_filtered/Bragg/7bragg_15dist.mat', 'C1')
    disp("DAS finished");
end
