% Contrast_Discrimination_Speed_Accuracy
% Preprocessing Steinemann dataset for unfold analyses


%% Specify file locations
% assumes that EEGlab is added to path already

% This is the directory of the folder for both data and m-files
prefix_ = ['~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Data/Steinemann_2018/'];
fileFolder = fullfile(prefix_, 'preprocessing_for_Unfold/analyses');
cd([fileFolder])

dataFolder = fullfile(prefix_, 'data');
% figData = dataFolder;
dataFolder_raw = fullfile(dataFolder, 'Raw data');

% added RF for saving out cleaned EEG files
dataFolder_out = [prefix_,'preprocessing_for_Unfold/data_clean'];

% Add Matlab packages to path
addpath(fileFolder)

% get subject foldernames to loop through for grabbing and concatenating
% data
tmpfoldnames = dir(dataFolder_raw);
foldnames = tmpfoldnames([tmpfoldnames(:).isdir]);
foldnames = foldnames(~ismember({foldnames(:).name}, {'.','..'}));

plotERPs=0;
processEEG =0;
%% Preprocessing
% This section rund on the full raw data set
% This data is not included in this repository and available upon request

fs=500; % sample rate

nchan = 93; % number of EEG channels
chanlocs = readlocs([fileFolder '/JBhead96_sym_EMG.loc']);    % electrode locations

headchans = [1:68 73:96];
headchanels = [1:68 73:97];

badblock = [];
badsub = [];
varsizes =[];
ref = [];
%% to do:
% loop through subjects and
allSubRTs = [];
for cursub= 1: length(foldnames)

    clear EEG curEEG

    subfold = foldnames(cursub).name; % get that partticipants folder with all the data.
    subpath = [dataFolder_raw,'/', subfold];

    fns = dir([subpath, '/*.vhdr']); % get all vhdr files
    fns = {fns(:).name}; % names only
    mfns = dir([subpath, '/*.mat']); % get all corresponding mat files
    ns = dir([subpath, '/*p*.mat']); % get files that are not behavioral files
    ns = {ns(:).name};
    mfns = {mfns(:).name}; % names only
    mfns = mfns(~ismember(mfns, ns)); % remove files that are not behavioral files
    B = regexp(fns,'\d*','Match'); % get block numbers to sort files
    b = cellfun(@(x)str2double(x),B); % convert into doubles for indexing
    bs = sortrows([b; 1:length(b)]', 1); % sort blocks with index columns
    fns = fns(bs(:,2)'); % now sort file names
    mfns = mfns(bs(:,2)');


    for curBlock = 1:length(fns)
        filename = fns{curBlock};
        matfilename = mfns{curBlock};

        % 1)  read in their data and merge all datasets
        curEEG = pop_loadbv(subpath,filename);    % EEGLAB function for opening BrainProducts EEG files

        curEEG.data = curEEG.data(headchans,:);% Isolate EEG data
        curEEG.data(nchan:end,:)=[];           % Empty out any channels beyond
        curEEG.data(end+1,:)=0;                % Add zeros for reference channel
        curEEG.nbchan = nchan;
        curEEG.chanlocs = chanlocs(1:nchan);
        curEEG.chanlocs(end).labels = 'FCz';


        load([subpath, '/', matfilename], 'resp')

        MksS = {'S  8', 'S136'};
        MksR = {'S 12', 'S 13'};
        Mks = {'S  8', 'S136', 'S 12', 'S 13'};
        % remove events that are not
        curEEG.event(~ismember({curEEG.event(:).type},  Mks))= [];
        % we expect that there is one event onset and one response for each
        % trial and that this matches the number of recorded response times


        curEventsS = {curEEG.event(ismember({curEEG.event(:).type},  MksS)).type};
        curEventsR = {curEEG.event(ismember({curEEG.event(:).type},  MksR)).type};
        rt = num2cell(repelem(resp.time,2))';
        accuracy = repelem(resp.response,2);

        [curEEG.event(:).rt] = deal(NaN);
        [curEEG.event(ismember({curEEG.event(:).type},  MksR)).code] = deal('Response');


        fprintf('Found %d stimulus onsets, %d responses and %d reaction times\n', length(curEventsS), length(curEventsR), length(rt))

        if length(rt)==(length(curEventsS)+ length(curEventsR))
            fprintf('Everything matches, adding RTs to event structure\n')
            [curEEG.event(:).rt] = rt{:};
            [curEEG.event(:).accuracy] = deal('correct');
            [curEEG.event(accuracy==3).accuracy] = deal('error');
        else
            fprintf('Oh, oh, something went wrong! Trying to fix...\n')
            if length(curEventsS)>length(curEventsR)
                % we identify instances when 2 stimulus triggers were sent
                % right after each other
                lns = find(ismember({curEEG.event(:).type},  MksS));
                repinds= lns(([lns(2:end),lns(end)+2]-lns)<2);
                % we then remove these instances
                curEEG.event(repinds)= [];
                curEventsS = {curEEG.event(ismember({curEEG.event(:).type},  MksS)).type};
                curEventsR = {curEEG.event(ismember({curEEG.event(:).type},  MksR)).type};
                fprintf('Found %d stimulus onsets, %d responses and %d reaction times\n', length(curEventsS), length(curEventsR), length(rt))

                if length(rt)==(length(curEventsS)+ length(curEventsR))
                    fprintf('Now everything matches, adding RTs to event structure\n')
                    [curEEG.event(:).rt] = rt{:};
                    [curEEG.event(:).accuracy] = deal('correct');
                    [curEEG.event(accuracy==3).accuracy] = deal('error');
                else

                    fprintf('Something is still off! Removing events.\n')
                    badblock = [badblock; curBlock];
                    badsub = [badsub; subfold];
                    varsizes = [varsizes; length(curEventsS), length(curEventsR), length(rt)];
                    curEEG.event(:) = [];

                end
            end

        end

        if curBlock > 1
            EEG = pop_mergeset(EEG, curEEG, 0); % merge the
        else
            EEG = curEEG;
        end

    end

    if processEEG
        % average reference
        EEG = pop_reref( EEG, []);

        % high pass filter at 40Hz (no low pass bc possibly slow potential)
        EEG = pop_eegfiltnew(EEG, 0.1, 40);

        pop_saveset(EEG, 'filename',sprintf('%s/%s_clean',dataFolder_out,subfold));
    end
    if plotERPs

        [sEEG, indices]=pop_epoch(EEG, MksS,[-0.2 1.8]);

        [rEEG, indices]=pop_epoch(EEG, MksR,[-1.8 0.2]);


        if cursub==1
            figure()
            set(gcf,'Color', [1 1 1])
        end
        subplot(2,length(foldnames),cursub)
        hold on
        title(sprintf('Participant %s', subfold))
        plot([-0.2:0.002:1.798], mean(sEEG.data( 53, :, :), 3))
        hline = refline([0 0]);
        set(hline,'Color','black','Linestyle',':','LineWidth',1);
        gridxy(0, 'Linestyle',':','LineWidth',1);
        xlabel('Time [s]')
        ylabel('s-locked Amplitude [µV]')
        subplot(2,length(foldnames),cursub+length(foldnames))
        plot([-1.8:0.002: 0.198], mean(rEEG.data( 53, :, :), 3))
        xlabel('Time [s]')
        ylabel('r-locked Amplitude [µV]')
        hline = refline([0 0]);
        set(hline,'Color','black','Linestyle',':','LineWidth',1);
        gridxy(0, 'Linestyle',':','LineWidth',1);
        set(gca,'TickDir','out', 'Box', 'off')
    end
    curSubRTs =  [EEG.event(ismember({EEG.event(:).type},  MksR)).rt];
    allSubRTs= [allSubRTs, curSubRTs];
end

save('AllSubRTs.mat', 'allSubRTs');
if plotERPs
    export_fig('AllSub_ERPs','-pdf','-painters');

end
