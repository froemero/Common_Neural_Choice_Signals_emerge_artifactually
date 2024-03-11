%% Specify file locations

% Location of Matlab packages
% matFolder = 'C:\Users\shadlenlab\Documents\MATLAB\';
% eeglab version: download eeglab (https://sccn.ucsd.edu/eeglab/download.php), code was originally run using eeglab14_1_1b
%eeglab_version = 'eeglab14_0_0b';
% Download toolbox norm_sphere_test for spericity testing: https://www.mathworks.com/matlabcentral/fileexchange/3694-sphertest
% Download bva-io toolbox for extracting data from Brain Vision Data Exchange format datasets from https://github.com/sccn/bva-io
% Download function CSD here: https://www.mathworks.com/matlabcentral/fileexchange/69399-current-source-density-csd

processEEG =0;

% This is the directory of the folder for both data and m-files
prefix_ = ['~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Data/Boldt_et_al_2019/'];
fileFolder = fullfile(prefix_, 'prep_Boldt_unfold/analyses');
cd([fileFolder])

dataFolder = fullfile(prefix_, 'squircle08_behav');
% figData = dataFolder;
dataFolder_raw = fullfile(prefix_, 'EEG');

% added RF for saving out cleaned EEG files
dataFolder_out = [prefix_,'prep_Boldt_unfold/data_clean'];

% Add Matlab packages to path
addpath(fileFolder)

allFiles = [dir(sprintf('%s/*.cnt', dataFolder_raw))];
allFiles = {allFiles(:).name}';
badsub = [];
varsizes =[];
for cursub = 1:length(allFiles)

    subID = allFiles{cursub}(11: end-12);

    % load that subject's behavioral data
    load([dataFolder,'/', 'squircle08_', subID, '.mat'], 'data');

    % remove any no-response trials (there's a bunch of NaN's in the RT
    % file and the number of entries doesn't match up with the number of
    % any triggers, but if we remove NaN's the RTs match the number of response triggers
    RTs = data.rt(~isnan(data.rt));
    length(data.rt(isnan(data.rt)));
    length(RTs);

    % load that subject's EEG data
    EEG = pop_loadcnt( [dataFolder_raw,'/',allFiles{cursub}], 'dataformat', 'auto', 'memmapfile', '');

    % remove EOG and M2
    EEG.data = EEG.data(1:32,:);
    EEG.chanlocs = EEG.chanlocs(1:32);
    EEG.nbchan = size(EEG.data,1);

    % match EEG and behavioral data - add RTs to event file
    MksS = [21, 22, 11, 12];
    MksR = [151, 152, 161, 162];
    Mks = [21, 22, 11, 12, 151, 152, 161, 162];

    % remove events that are not stimulus or response triggers
    EEG.event(~ismember([EEG.event(:).type],  Mks))= [];
    length(EEG.event);

    % we expect that there is one event onset and one response for each
    % trial and that this matches the number of recorded response times
    curEventsS = [EEG.event(ismember([EEG.event(:).type],  MksS)).type];
    curEventsR = [EEG.event(ismember([EEG.event(:).type],  MksR)).type];

    fprintf('Found %d stimulus onsets, %d responses and %d reaction times\n', length(curEventsS), length(curEventsR), length(RTs))


    % there are a bunch of duplicate stimulus onsets, so we look for where
    % multiple stimulus onsets occur
    lns = find(ismember([EEG.event(:).type],  MksS));
    repinds= lns(([lns(2:end),lns(end)+2]-lns)<2);
    % we then remove these instances
    EEG.event(repinds)= [];
    curEventsS = [EEG.event(ismember([EEG.event(:).type],  MksS)).type];
    fprintf('Found %d stimulus onsets, %d responses and %d reaction times\n', length(curEventsS), length(curEventsR), length(RTs))

    if (length(curEventsS)==length(curEventsR)) && (length(curEventsR)==length(RTs))

        rt = num2cell(repelem(RTs,2))';
        [EEG.event(:).rt] = deal(NaN);
        [EEG.event(:).accuracy] = deal('correct');

        [EEG.event(ismember([EEG.event(:).type],  MksR)).code] = deal('Response');
        [EEG.event(ismember([EEG.event(:).type],  [161, 162])).accuracy] = deal('error');
        [EEG.event(find(ismember([EEG.event(:).type],  [161, 162]))-1).accuracy] = deal('error');
        [EEG.event(ismember([EEG.event(:).type],  MksS)).code] = deal('Stimulus');
        [EEG.event(:).rt] = rt{:};

        curSubRTs =  [EEG.event(ismember([EEG.event(:).type],  MksR)).rt];
        allSubRTs= [allSubRTs, curSubRTs];

        if processEEG ==1
        % high pass filter at 40Hz (no low pass bc possibly slow potential)
        EEG = pop_eegfiltnew(EEG, [], 40);
        pop_saveset(EEG, 'filename',sprintf('%s/%s_clean',dataFolder_out,subID));
        end
    else
        fprintf('Oh, oh, something went wrong!\n')
        badsub = [badsub; str2double(subID)];
        varsizes = [varsizes; length(curEventsS), length(curEventsR), length(RTs)];
    end



end
save('AllSubRTs.mat', 'allSubRTs');
