%% % goal: generate table of format like the files Bene sent:

% subject, type (epoched, deconv), basisname (Stimulus, Response),
% condition (here we have 4, not 2 like in the other datasets), time, yhat


% this will be 40*2*4*timepoints(this is different for s & r, so let's
% maybe not preallocate all of it.

%% prep: get all relevant things
PATH = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Analyses/Matlab/UCAP/';

%% get channel locations

load(strcat(PATH, 'Second_level/Export/chanlocs.mat'));
addpath(strcat('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/', 'helpfunctions'))
addpath(genpath('~/Dropbox (Brown)/CLPS-ShenhavLab/Resources/EEG_resources/EEGfunctions/'));

for n = 1:size(chanlocs,2)

    expression = [chanlocs(n).labels '=' sprintf('%d',n) ';'];
    eval(expression)
end

%%

% load behavioral Data
load(strcat(PATH, 'Second_level/Export/behavTable.mat'));

behavTable.Quality = zeros(size(behavTable.VPNummer));
behavTable.Quality(strcmp(behavTable.n_b,'normal')==1)=1;
behavTable.DevP = zeros(size(behavTable.VPNummer));
behavTable.DevP(strcmp(behavTable.DeviantPosRL,'re')==1)=1;
behavTablecor = behavTable(behavTable.accuracy==1,:);
%%

% load regcoeffs and ROI data for that analysis
load('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Analyses/Matlab/UCAP/Second_level/Export/regcoeffs.mat')
load('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Analyses/Matlab/UCAP/Second_level/Export/DeviantP_Quality_DPbyQuality_0.005/ROI_data.mat')

%% subsetting only Pz
% (electrodes, time-points, condition, subject)
sbasicCoeffs.data= sbasicCoeffs.data(Pz, [find(sbasicCoeffs.times== -0.2) : find(sbasicCoeffs.times== 0.998)],:,:);
rbasicCoeffs.data= rbasicCoeffs.data(Pz, [find(rbasicCoeffs.times== -1) : find(rbasicCoeffs.times== 0.2)],:,:);

%%
[~,idx] = unique(behavTable.VPNummer,'first');
vps = behavTable.VPNummer(sort(idx));



%%
clear yHat
clear allYHats
allYHats = [];
lastn = 0;
curbetas = sbasicCoeffs.data; % let's do one first - first Valuation component
for subj = 1:length(vps)
    s_id = vps(subj);
    clear yHat
    bSubData = behavTable(behavTable.VPNummer==s_id & behavTable.accuracy~= 0,:); % get a subjects data, I think I already removed non-responses, but whatever

    xMat=meanCentX([ones(size(bSubData,1),1),bSubData.DevP , bSubData.Quality, bSubData.DevP.*bSubData.Quality]); %zeros(size(bSubData,1),1) now I'm setting the intercept and the variable of interest to zero

    for channel = 1: size(curbetas,1)
        for t = 1: size(curbetas,2)
            beta = squeeze(curbetas(channel,t,:, subj));
            yHat(channel, t, :)= [xMat]*beta;
        end
    end

    allYHats(:,:,lastn+1:lastn+size(yHat,3)) = yHat;
    lastn = size(allYHats,3);
end
Residuals = allYHats;





