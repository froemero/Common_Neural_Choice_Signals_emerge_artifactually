%% Attention, stuff has been altered for no filter analyses

%% Info for segmentation analyses:

%  time-lock to stimulus onset

PATH = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Data/';
addpath('~/Dropbox (Brown)/CLPS-ShenhavLab/Resources/EEG_resources/EEGfunctions/')

doSeg=1;

Trigs =[10,20];
Mks = cellstr(num2str(Trigs(:)));
Mks=strcat({'S'}, {' '}, Mks); % adding S to Markers for segmentation

%% for response-locked segementation:
load(sprintf('%sExport/allSubDataTable.mat', PATH))

%% approach: Get all condition triggers
% just NAN the artifact trials and paste all subject files together.
%% run only when needed, adds EEGLab

% addpath('N:/Software/eeglab13_5_4b')
%% preparation
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%%
SOURCEFILES = dir(strcat(PATH, 'EEG_clean/*_cleanNF.set')); %all cleaned files
SUBJECTS = 1:numel(SOURCEFILES);
%% Combine pre- and post response segmentation to one

DAT=nan(65,2100,40*120); % matrix: electrode * time points * (participants*trials)
DATR=nan(65,1400,40*120); % matrix: electrode * time points * (participants*trials)
PR=0;
%DATPR=nan(65,500,40*120); % for ERN/Pe -200: 800 post response
contmat=zeros(40,3, length(Mks)); % check number of trials found per participant (before & after AR)

SubAllLatencies= [];
%%
for s =1:numel(SOURCEFILES)
    %%
    
    
    vpn=SOURCEFILES(s).name(1:4); % read out of dataset name
    s_id = str2double(SOURCEFILES(s).name(3:4)); % reads in subject ID and converts it to double (does not work for '01')
    fprintf('processing participant number %d\n', s_id)
    contmat(s,1,:)=s_id;
    pcnt=(s*120)-119; % counter line ntrials; bad luck with presentation programming: some triggers are repeated with different meaning;
    minRT = min(allSubDataTable.RTeval1(allSubDataTable.SubNum==str2double(vpn)));
    %%
    if ~isempty(minRT)
        
        EEG = pop_loadset('filename',sprintf('%s_clean.set',vpn),'filepath',sprintf('%sEEG_clean/',PATH));
        
        EEGO=EEG;
       
        %%

        vpdouble = str2double(vpn);
        % if response is nan, add 4000/2
        subRTs = allSubDataTable.RTeval1(allSubDataTable.SubNum==vpdouble)*1000;
        subLatencies= nan(length(subRTs),1);
        subLatenciesCHECK= nan(length(subRTs),1);
        subRTs(isnan(subRTs)) = 4000;
        lcnt=1;
        fprintf('%d\n', size(EEG.event,2))
        for eventnumber=11:2:size(EEG.event,2) % the first two triggers are new segment and acticap, the next 8 are practice trials, we skip that 
         
         if s_id==39 % subject 39 has 1 trigger less than the others because that dillhole started before they should. So here I need to subtract 1 from each; but that subject has the block-start and end trigger... sooo...
            if lcnt <= length(subRTs)
           subLatencies(lcnt) = EEG.event(eventnumber+1).latency-EEG.event(eventnumber).latency;  
           EEG.event(eventnumber+1).latency =  EEG.event(eventnumber).latency + round(subRTs(lcnt)/2);
           subLatenciesCHECK(lcnt)= EEG.event(eventnumber+1).latency-EEG.event(eventnumber).latency;
            end
         elseif size(EEG.event,2)>250
             if lcnt <= length(subRTs)
           subLatencies(lcnt) = EEG.event(eventnumber+2).latency-EEG.event(eventnumber+1).latency;  
           EEG.event(eventnumber+2).latency =  EEG.event(eventnumber+1).latency + round(subRTs(lcnt)/2);
           subLatenciesCHECK(lcnt)= EEG.event(eventnumber+2).latency-EEG.event(eventnumber+1).latency;
             end
         else
             % subject 5012 has 252 triggers... one is S60... I don't send
             % that...
           subLatencies(lcnt) = EEG.event(eventnumber+1).latency-EEG.event(eventnumber).latency;  
           EEG.event(eventnumber+1).latency =  EEG.event(eventnumber).latency + round(subRTs(lcnt)/2);
           subLatenciesCHECK(lcnt)= EEG.event(eventnumber+1).latency-EEG.event(eventnumber).latency;
         end
          lcnt=lcnt+1;  
        end 
        CHECK= [subLatencies, subLatenciesCHECK];
        CHECK2 = CHECK*2;
        CHECK2(:,3) = CHECK2(:,1)- CHECK2(:,2);
        SubAllLatencies=[SubAllLatencies;CHECK2];
        %%
         EEGC=EEG; % save copy in temp memory
         if doSeg== 1
        %% get stim-locked data
        for mnr= 1: length(Mks)
            if mnr ==1
                fprintf('s-locked segmentation of participant number %d\n', s_id)
                IL =[-0.2 4]; % 200 ms pre- stimulus + entire trial length
                DP = size(DAT,2); % number of data points in segment
                [EEG, indices]=pop_epoch(EEG, Mks(mnr),IL);
                EEG = pop_rmbase( EEG, [-200 0]); % Baseline correction relative to stimulus onset
            else
                fprintf('r-locked segmentation of participant number %d\n', s_id)
                if PR==1
                IL =[-0.2 0.8]; % max RT - 300 ms, so we include this rather s-locked early P3 for long RT-files, too
                DP = size(DATPR,2);   
                [EEG, indices]=pop_epoch(EEG, Mks(1),[(-4.3+minRT) 5.3]);
                else
                %IL =[-3.7 0.2]; % max RT - 300 ms, so we include this rather s-locked early P3 for long RT-files, too
                IL =[-2 0.8]; % for both pre and post response
                DP = size(DATR,2);
                % in order to have a reasonable baseline, we first do a regular
                % stimulus locked segmentation with a large time window, so we can fit in all possible
                % response durations with the same segmentation later on. For
                % each individual we determine the maximum window size by
                % computing the smallest RT and adding it to -4 (plus 100 ms to
                % prevent rounding issues):
                % stimlocked segmentation with variable interval based on participant performance
                %[EEG, indices]=pop_epoch(EEG, Mks(1),[(-4.3+minRT) 4.5]); % much bigger window required than expected... response triggers and/or triggers in general seem not very accurate... 
                [EEG, indices]=pop_epoch(EEG, Mks(1),[(-2.5+minRT) 5.5]); % much bigger window required than expected... response triggers and/or triggers in general seem not very accurate... 

                end
                EEG = pop_rmbase( EEG, [-200 0]); % Baseline correction relative to stimulus onset
                % sub epoch requires only the marker position relativ to the
                % current 0 event (here 1) and a time interval (IL).
                [EEG, indices]=sub_epoch(EEG, 1, IL);
                %[EEG, indices]=pop_epoch(EEG, Mks(mnr),IL);
            end
            % select all trials based on object onset markers
            % no further segmentation required --> matching with logfiles
            contmat(s,2,mnr)=EEG.trials;
            %%
            % AR: Ansatz: Trials werden nicht gel?scht, sondern nur die
            % unzul?ssigen identifiziert und dann sp?ter nicht mit ausgelesen.
            % (Irej: index of rejected trials)
            
            [EEGn, Irej] = pop_eegthresh(EEG,1,[1:60 65] ,-150,150, IL(1),IL(2)-0.002,1,0); % artifact rejection based on all but the ocular electrodes
            
            %% %   >> [rej rejE] = rejtrend( signal, winsize, maxslope, minR, step);
            
            [rej, rejE] = rejtrend(EEG.data,DP,50,0.3,1); % worked so far ony all data, but check again when you have more subjects (1 trial excluded for 5000)
            % (rej: vector with 0 & ones... bad!)
            
            frej= find(rej==1);
            Out=unique(sort([Irej,frej]));
            if isempty(Out)
                remout=0; 
            else
            remout = length(Out);
            end
            fprintf('%d trials removed alltogether\n', remout);
            %%
            EEG.data(:,:,Out) = nan; % invalidate data with artifacts
            %%
            % define which trials are not part of the actual task but
            % training --> ptrials is number of practice trials + 1 (the
            % first trial after practice
            if s_id == 39
                ptrials = 4; % debug for subject who started task before recording was on, so that one practice trial is missing
            else
                ptrials = 5;
            end
                
            if mnr ==1
                %EEG = pop_rmbase( EEG, [-200 0]); % Baseline correction
                DAT(:,:,pcnt:pcnt+size(EEG.data,3)-ptrials)=EEG.data(:,:,ptrials:size(EEG.data,3));% excludes practice trials
            else
                if PR==1
                DATPR(:,:,pcnt:pcnt+size(EEG.data,3)-ptrials)=EEG.data(:,:,ptrials:size(EEG.data,3));% excludes practice trials
                else
                DATR(:,:,pcnt:pcnt+size(EEG.data,3)-ptrials)=EEG.data(:,:,ptrials:size(EEG.data,3));% excludes practice trials
                end
                
            end
            contmat(s,3,mnr)=(EEG.trials-size(Out,2));
            EEG=EEGC; % get original EEG back to do the second segmentation
        end
        end
    else
        fprintf('run analyze BASB_copy first!\n')
    end
end

%%
fprintf('%d participants in dataset... saving data\n', s)
save(sprintf('%sExport/DATnew.mat', PATH), 'DAT', '-v7.3');
save(sprintf('%sExport/DATRnew.mat', PATH), 'DATR', '-v7.3');
%save(sprintf('%sExport/DATPR.mat', PATH), 'DATPR', '-v7.3');