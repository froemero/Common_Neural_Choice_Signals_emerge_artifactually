% Export time window average or time series to access in R for analyses  


PATH = '/Users/Romy/Dropbox (Brown)/ShenhavLab/EEG_ressources/Experiments/BASB_EEG/'; %edit to your needs
EXPORTPATH ='Data/Export/'; % Export folder (segmented data will be saved there)
srate=500; % sampling rate
bl=200; %length of prestimulus interval in ms

%% 1) load EEG data (if necessary)

load(strcat(PATH,EXPORTPATH,'DAT2new.mat'), 'DAT2'); % load whatever your data is
load(strcat(PATH,EXPORTPATH,'DATR2new.mat'), 'DATR2');


%% EXPORT CLUSTERS FOR LMM ANALYSES

% choice component r-locked: tw: 990 - 500 --> 
BFS = DATR2; % set BFS to whatever your data is
nchans = size(BFS,1);

%% export mean cluster time
% clumsy way to deal with negative time:
% use time until zero and then subtract actual time
%Cluster time: -954 ms : -554 ms 
mint=2000-954; %export start in ms
maxt=2000-554; %export end in ms
ERPNAME='ChoiceClust'; % name of the component you want to export
bl=0;
% the rest is automatic
EXPERP=nanmean(BFS(:,(mint+bl)*srate/1000:(maxt+bl)*srate/1000,:),2); % adapt for multiple sampling rates

EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');


%% export Valuation component
% choice component r-locked: tw: 990 - 500 --> 
BFS = DAT2; % set BFS to whatever your data is
nchans = size(BFS,1);

%% export mean cluster time
% clumsy way to deal with negative time:
% use time until zero and then subtract actual time
%Cluster time: 702 ms : 782 ms 
mint=702; %export start in ms
maxt=782; %export end in ms
ERPNAME='ValClust'; % name of the component you want to export
bl=200;
% the rest is automatic
EXPERP=nanmean(BFS(:,(mint+bl)*srate/1000:(maxt+bl)*srate/1000,:),2); % adapt for multiple sampling rates

EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');


%% for frontal cluster 
% clumsy way to deal with negative time:
% use time until zero and then subtract actual time
mint=686; %export start in ms
maxt=782; %export end in ms
ERPNAME='ValClustf'; % name of the component you want to export
bl=200;
% the rest is automatic
EXPERP=nanmean(BFS(:,(mint+bl)*srate/1000:(maxt+bl)*srate/1000,:),2); % adapt for multiple sampling rates

EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');


%% try the CPP the Pisauro way (average -700 to -200)
% here I exported the naned out non trial data
% zero is 1000 (2000ms);  -700 is (2000 - 700)/2 (650), -200 is
% (2000-200)/2 (900)

mint=-700; %export start in ms
maxt=-200; %export end in ms
ERPNAME='CPPPisauro'; % name of the component you want to export
bl=2000;
% the rest is automatic
EXPERP=nanmean(BFS(:,(mint+bl)*srate/1000:(maxt+bl)*srate/1000,:),2); % adapt for multiple sampling rates

EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',abs(mint),abs(maxt)) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,abs(mint),abs(maxt)), sprintf('%s%d%d',ERPNAME,abs(mint),abs(maxt)), '-v6');

