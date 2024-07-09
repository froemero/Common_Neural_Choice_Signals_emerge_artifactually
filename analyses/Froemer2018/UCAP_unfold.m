
%%

PATH = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Data/';
addpath('~/Dropbox (Brown)/CLPS-ShenhavLab/Resources/EEG_resources/EEGfunctions/')
addpath(strcat('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/', 'helpfunctions'))
        
%% load relevant data and merge PCs with rest of behavior
% load behavior
load(sprintf('%sUCAP_Export/behav.mat', PATH))
behavTable= struct2table(a1);
       % save( strcat(PATH, 'Second_level/Export/behavTable.mat'), 'behavTable')
%% start unfold
addpath(genpath('~/Dropbox (Brown)/CLPS-ShenhavLab/Resources/EEG_resources/unfold/'))
init_unfold
eeglab;
%% Goal: run unfold analyses on all participants 

%%
vps =unique(behavTable.VPNummer);

% SOURCEFILES = dir(strcat(PATH, 'EEG_clean/*_cleanNF.set')); %all cleaned files
for Subject = 1:length(vps)

    s_id = vps(Subject); % reads in subject ID and converts it to double (does not work for '01')
    if s_id <10
     vpn=strcat('0',num2str(s_id));   % add zeros for single digigts so we can actually read in the files
    else
    vpn=num2str(s_id); % read out of dataset name
    end
    fprintf('processing participant number %d out of 40\n', Subject)

%% load one dataset
        EEG = pop_loadset('filename',sprintf('%s_clean.set',vpn),'filepath',sprintf('%sUCAP_clean/',PATH));
        %% recode markers to make my life easier
        
        [EEG.event(find(strcmp({EEG.event.type}, 'S 10')==1)).type] = deal('S 500'); % getting preexisting 'S 10' out of the way 
        
        % try resetting all stimulus onsets to S10
        %% these are all stimulus onset triggers
% 

Trigs =[201:208, 211:218]; % specify stimulus onset trigger(s)
Mks = cellstr(num2str(Trigs(:)));
Mks=strcat('S', Mks); % adding S to Markers for segmentation;

for currtrigg =  1: length(Mks)
    currMk = Mks{currtrigg};
[EEG.event(find((strcmp({EEG.event.type}, currMk)==1))).type] = deal('S 10');

end


 %% ok, so we have the correct number of stim onsets
allStims = find(strcmp({EEG.event.type}, 'S 10')==1);
isStim = (strcmp({EEG.event.type}, 'S 10')==1); 
isResp = (strcmp({EEG.event.type}, 'S 99')==1 | strcmp({EEG.event.type}, 'S 98')==1 );
isStimShift = zeros(size(isStim));
isStimShift(2:end)= isStim(1:end-1); 
% true responses, where previous event was a stim onset
RespIndex = find(isResp==1 & isStimShift==1);
% replace only responses after stimuli with
% the marker we want

[EEG.event(RespIndex).type] = deal('S 20');
%% get event info from table and add to event structure in EEG
        % add continuous predictors to EEG.event structure
        SubData=behavTable(behavTable.VPNummer==s_id,:);
% we need:
% 1) blurr vs intact (col 5)
% 2) DeviantPosRL (col 12)
% 3) err (col 32)

% % scale values where required and mean center regressors
%                 % THIS IS FUCKING AWESOME!
%         SubData{:,[7:13 18 19]}= SubData{:,[7:13 18 19]}./10;
%         SubData{:,[7:13 18 19 44:45]}= meanCentX(SubData{:,[7:13 18 19 44:45]});

        ofieldNames = fieldnames(EEG.event);
        % this is the info we want to add to the EEG 44 and 45 are the PCs
        colindx = [5 12 32];
        
%         if Subject < 6 
%             discarter = 10;
%         else
%             discarter = 11;
%         end
%%         
        for IVs = 1: length(SubData(:,colindx).Properties.VariableNames)
                    cnt=0;

            currname= SubData(:,colindx).Properties.VariableNames{IVs};
            for iii= 1:length(EEG.event)
                EEG.event(1).(currname)=0;
                if strcmp(EEG.event(iii).type, 'S 10')
%                     if iii<discarter
%                         EEG.event(iii).(currname)=0;
%                     else
                        cnt=cnt+1;
                        EEG.event(iii).(currname)=SubData{cnt,colindx(IVs)};
                        if strcmp(EEG.event(iii+1).type, 'S 20')
                            EEG.event(iii+1).(currname)=SubData{cnt,colindx(IVs)};
                        end
%                     end
                    
                elseif (isempty(EEG.event(iii).(currname)))
                    EEG.event(iii).(currname)=0;
                    
                end
                
            end
        end
        %% exclude events that are not part of the actual data
        % we don't have a lot of errors here, so we can't estimate them
        % properly and hence exclude them
[EEG.event(find((strcmp([EEG.event(:).err], 'wrong')==1))).type] = deal('S 00');

%%
EEG = eeg_checkset(EEG);
cfgDesign = [];
cfgDesign.eventtypes = {'S 10','S 20'};
%%
%cfgDesign.formula = {'y ~ 1+AvBid+ MaxvMinBid','y~1+AvBid+ MaxvMinBid'};
cfgDesign.formula = {'y ~ 1+ cat(n_b) * cat(DeviantPosRL) ','y~1+ cat(n_b) * cat(DeviantPosRL)'};
EEG = uf_designmat(EEG,cfgDesign);

%%
cfgTimeexpand = [];
cfgTimeexpand.timelimits = [-1,1];
%cfgTimexpand.method = 'stick';
cfgTimexpand.method = 'fourier'; %% hat eigentlich stick gemacht!!!! 
cfgTimexpand.timeshiftparam = 20; 
EEG = uf_timeexpandDesignmat(EEG,cfgTimeexpand);

% plot design matrix if desired
% uf_plotDesignmat(EEG)
% uf_plotDesignmat(EEG,'timeexpand',1)

%% Artifact detection
winrej = uf_continuousArtifactDetect(EEG,'amplitudeThreshold',250);

%% Artifact rejection
EEG = uf_continuousArtifactExclude(EEG,struct('winrej',winrej));

%% Fit deconvolution
EEG= uf_glmfit(EEG);

%% epoching
EEG = uf_epoch(EEG,struct('winrej',winrej,'timelimits',cfgTimeexpand.timelimits));

%% fit massunivariate on epoched data without deconvolution
EEG = uf_glmfit_nodc(EEG); 
%% make results file
ufresult= uf_condense(EEG);

%% Save out results for each subject 
save(sprintf('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Analyses/Matlab/UCAP/UCAP_UF_results/%s_ufresult.mat', vpn), 'ufresult', '-v7.3');
%display(ufresult)

% %%
% uf_plotParam(ufresult,'channel',21);
% %%
% g = uf_plotParam(ufresult,'channel',50,'deconv',1,'baseline',[ufresult.times(1) 0]);
% 
% g = uf_plotParam(ufresult,'channel',50,'deconv',0,'baseline',[ufresult.times(1) 0]);
% 
% uf_plotParam(ufresult,'channel',21);
% 
% g = uf_plotParam(ufresult,'channel',21,'deconv',1);

end

