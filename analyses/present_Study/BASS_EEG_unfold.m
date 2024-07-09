
%%

PATH = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASS_EEG/Data/';
addpath('~/Dropbox (Brown)/CLPS-ShenhavLab/Resources/EEG_resources/EEGfunctions/')
addpath(strcat('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/', 'helpfunctions'))
        
%% load relevant data and merge PCs with rest of behavior
% load behavior
load(sprintf('%sexport/allSubDataTable505.mat', PATH))
% load factor score table
load(sprintf('%sexport/PCAfactorScoreTable.mat', PATH))

% remove subjects without EEG
allSubDataTable = allSubDataTable((allSubDataTable.SubNum ~= 5000 & allSubDataTable.SubNum ~= 5022), :);

% create merged data frame
allSubDataTablePC = allSubDataTable;
allSubDataTablePC.AppraisalPC(allSubDataTable.Choice~=-1,:) = scoretable.Valuation;
allSubDataTablePC.ChoicePC(allSubDataTable.Choice~=-1,:) = scoretable.Choice;


%% start unfold
addpath(genpath('~/Dropbox (Brown)/CLPS-ShenhavLab/Resources/EEG_resources/unfold/'))
init_unfold
eeglab;
%% Goal: run unfold analyses on all participants that have Subeval data. Save out the results and run CBPT a la Matt on the resulting betas


%%
vps = unique(allSubDataTablePC.SubNum);

% SOURCEFILES = dir(strcat(PATH, 'EEG_clean/*_cleanNF.set')); %all cleaned files
for Subject = 1:length(vps)

    s_id = vps(Subject); % reads in subject ID and converts it to double (does not work for '01')
    vpn=num2str(s_id); % read out of dataset name
    fprintf('processing participant number %d\n', s_id)

    %% load one dataset
    EEG = pop_loadset('filename',sprintf('BASS_%s_clean.set',vpn),'filepath',sprintf('%sclean/',PATH));
    %% get event info from table and add to event structure in EEG
    % add continuous predictors to EEG.event structure
    SubData=allSubDataTablePC(allSubDataTablePC.SubNum==s_id,:);

    % scale values where required and mean center regressors
    % THIS IS FUCKING AWESOME!
    %SubData{:,[7:13 18 19]}= SubData{:,[7:13 18 19]}./10; % we don't
    %actually need any of the other variables for our current purposes
    %atm.
    %SubData{:,[7:13 18 19 44:45]}= meanCentX(SubData{:,[7:13 18 19 44:45]});
    SubData{:,[size(SubData, 2)-1 size(SubData, 2)]}= meanCentX(SubData{:,[size(SubData, 2)-1 size(SubData, 2)]});

    ofieldNames = fieldnames(EEG.event);
    % this is the info we want to add to the EEG 44 and 45 are the PCs
    %colindx = [7:13 18 19 44:45];
    colindx = [3 size(SubData, 2)-1 size(SubData, 2)]; % adding trialnumber for control purposes


    %%
    for IVs = 1: length(SubData(:,[3 size(SubData, 2)-1 size(SubData, 2)]).Properties.VariableNames)
        cnt=0;

        currname= SubData(:,[3 size(SubData, 2)-1 size(SubData, 2)]).Properties.VariableNames{IVs};
        for iii= 1:length(EEG.event)
            EEG.event(1).(currname)=0;
            if strcmp(EEG.event(iii).type, 'S  1')
                if cnt==0
                discarter = iii-1;
                end
                cnt=cnt+1;
                %                     if iii<discarter
                %                         EEG.event(iii).(currname)=0;
                %                     else
                EEG.event(iii).(currname)=SubData{cnt,colindx(IVs)};
            elseif strcmp(EEG.event(iii).type, 'S 10')
                EEG.event(iii).(currname)=SubData{cnt,colindx(IVs)};
            elseif strcmp(EEG.event(iii).type, 'S 20')
                EEG.event(iii).(currname)=SubData{cnt,colindx(IVs)};
                
            elseif strcmp(EEG.event(iii).type, 'S 99')
                EEG.event(iii).(currname)=SubData{cnt,colindx(IVs)};
                %                     end
            elseif (isempty(EEG.event(iii).(currname)))
                EEG.event(iii).(currname)=0;

            end

        end
    end
    %% exclude events that are not part of the actual data
    EEG.event = EEG.event(discarter+1:end); % this might need to be set up flexibly due to 1 participant who started early and misses practice trials in recording
    % for some reason that still matched for the participant... should
    % check other subject.
    %%
    EEG = eeg_checkset(EEG);
    cfgDesign = [];
    %cfgDesign.eventtypes = {'S  1', 'S 10','S 20'};
    cfgDesign.eventtypes = {'S  1', 'S 20'};
    %%
    %cfgDesign.formula = {'y ~ 1+AvBid+ MaxvMinBid','y~1+AvBid+ MaxvMinBid'};
    %cfgDesign.formula = {'y ~ 1+AppraisalPC+ ChoicePC','y ~ 1','y~1+AppraisalPC+ ChoicePC'};
    cfgDesign.formula = {'y ~ 1+AppraisalPC+ ChoicePC','y~1+AppraisalPC+ ChoicePC'};

    EEG = uf_designmat(EEG,cfgDesign);


    %%
    cfgTimeexpand = [];
    cfgTimeexpand.timelimits = [-2,2];
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
    save(sprintf('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Analyses/Matlab/Unfold_Analyses_BASS/UF_Results/%s_ufresult.mat', vpn), 'ufresult', '-v7.3');
    %display(ufresult)

    % %%
    % uf_plotParam(ufresult,'channel',50);
    % %%
    % g = uf_plotParam(ufresult,'channel',50,'deconv',1,'baseline',[ufresult.times(1) 0]);
    %
    % g = uf_plotParam(ufresult,'channel',50,'deconv',0,'baseline',[ufresult.times(1) 0]);
    %
    % uf_plotParam(ufresult,'channel',21);
    %
    % g = uf_plotParam(ufresult,'channel',21,'deconv',1);

end

