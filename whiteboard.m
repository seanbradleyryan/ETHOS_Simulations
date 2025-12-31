clear all; 
wd = '/mnt/weka/home/80030361/ETHOS_Simulations/1194203/Pancreas/Session_1/sct';
wd = '/home/80030361/ETHOS_Simulations/EthosExports/1194203/Pancreas/Session_1/sct'
rtplan = fullfile(wd,'RTPLAN.1.2.246.352.800.5261851517227480225.9292172278832325515'); 

rtplan = dicomread(rtplan); 