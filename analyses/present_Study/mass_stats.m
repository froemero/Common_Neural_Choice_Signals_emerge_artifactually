%% To do mass-univariate ala Matt (Nassar)
% script adapted from Matt's cannon task script
% if anything breaks, contact romy.froemer@googlemail.com

% load data
PATH = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/';
EXPORTPATH = 'Data/Export/';
figDir= strcat(PATH, 'Docu/Figures/');
load(sprintf('%sData/Export/allSubDataTable.mat', PATH))

load(sprintf('%sData/Export/DATnew.mat', PATH))
load(sprintf('%sData/Export/DATRnew.mat', PATH))
% see bottom of this script for how these are generated
load(sprintf('%sData/Export/DAT2new.mat', PATH))
load(sprintf('%sData/Export/DATR2new.mat', PATH))
% load factor score table; making code below
load(sprintf('%sData/Export/PCAfactorScoreTable.mat', PATH))
% load this code to make it is at the end of the script
load(strcat(PATH, 'Data/Export/connectionMat.mat'))
noEOG=1; % with or without EOG channels
doViewing = 0; % if you want to include viewing location as a time-varying regressor.
% note that the script doesn't currently do that, though...

if doViewing % if I happen to want to do anything with the EOG data, I can use this
    load(sprintf('%sData/Export/RnSFixationLoc.mat', PATH))
end
%% loop through subjects
%% get unique subject names for later processing and preselection of partially existing data
[~,idx] = unique(allSubDataTable.SubNum,'first');
vpss = allSubDataTable.SubNum(sort(idx));
vps=vpss(vpss~=5028);% excluding subject 5028 because of missing Liking data & vpss~=5024

%% get channel locations and channel numbers (for plotting)
nchans = size(DAT,1);
load(sprintf('%sData/Export/chanlocs.mat', PATH));
for n = 1:nchans

    expression = [chanlocs(n).labels '=' sprintf('%d',n) ';'];
    eval(expression)
end
%% in this script we don't use EOG, so we exclude those channels
nchans=61; % reset number of channels
chanlocs= chanlocs([1:60,65]); % adjust chanlocs
Cz=61; % adjust channumber for Cz

%% compute value left - value right

allSubDataTable.Vlminr = allSubDataTable.Valuel - allSubDataTable.Valuer;

%% compute time to deadline
allSubDataTable.tToDL = 4-allSubDataTable.RTeval1;

%% generate path to helper functions
addpath(strcat(PATH, 'helpfunctions'))

% contains:
% meanCentx - for centering predictors
% nanzscore - the name says it
% zScoreX   - I assume that's being used by some of the other functions?
% getEEG_clusterSize - to obtain positive and negative clusters during
% permutation test

%%
% set parameters

corrType='spearman';
whichdata = "s"; % "r" or "s"
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
doSigTesting=1;
weightedT=1;
byHand=0;
numPerm=1000;
ocbColors=[0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167]./256;
cbColors= [ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors];
inc=2; % right now we have a sampling rate of 500Hz, if I take every 4th of this we sample every 8 ms
% so I would have a sampling rate of 125 (that's pretty low... on the other
% hand I'm only looking at slow changes... might be ok...
%% We need this to match the factor scores (which don't have subject 5028 or missing response trials)
suballSubDataTable= allSubDataTable(allSubDataTable.SubNum~=5028 & allSubDataTable.Choice1~=-1,:);
%% Here we loop through a vector that spacifies which analyses should be run

ModsToRun = [1:4];
AllThreshs= [.001, .005]; % cluster forming threshold. We confirm that results hold at more conservative threshold as well 
for grabMod = 1:length(ModsToRun) % loop through models I want to run

    for currthresh= 1:length(AllThreshs) % loop through thresholds I want to test
        thresh=AllThreshs(currthresh);
        segs=['r';'s';'p']; % add pr?
        clear 'sbasicCoeffs'
        clear 'rbasicCoeffs'
        clear 'prbasicCoeffs'

        tw=0;
        twl='TW2s';
        isROI=0;
        isEOG =0;
        ROIelecs = [CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4];
        ROIname = 'absEOG';


        if ModsToRun(grabMod) == 1
            % 1) SUBEVAL MODEL
            sbasicCoeffs.labels={'intercept', 'Liking', 'Anxiety', 'Confidence'};
            rbasicCoeffs.labels={'intercept', 'Liking', 'Anxiety',  'Confidence'};
            prbasicCoeffs.labels={'intercept', 'Liking', 'Anxiety',  'Confidence'};
            usefactorscores=0;
        elseif  ModsToRun(grabMod) == 2
            % 2) CHOSEN UNCHOSEN MODEL
            sbasicCoeffs.labels={'intercept', 'Chosen', 'Unchosen' }; %
            rbasicCoeffs.labels={'intercept', 'Chosen', 'Unchosen'}; %
            prbasicCoeffs.labels={'intercept', 'Chosen', 'Unchosen'};
            usefactorscores=0;
        elseif  ModsToRun(grabMod) == 3
            % 3) OV RV MODEL
            sbasicCoeffs.labels={'intercept', 'VD', 'OV'};
            rbasicCoeffs.labels={'intercept', 'VD', 'OV'};
            prbasicCoeffs.labels={'intercept','VD', 'OV'};
            usefactorscores=0;
        elseif  ModsToRun(grabMod) == 4
            % 4) PCA
            sbasicCoeffs.labels={'intercept', 'PCValuation', 'PCChoiceDifficulty'};
            rbasicCoeffs.labels={'intercept', 'PCValuation', 'PCChoiceDifficulty'};
            prbasicCoeffs.labels={'intercept','PCValuation', 'PCChoiceDifficulty'};
            usefactorscores=1;

        end
        SsT = 'G1500';


        if weightedT
            out_dir = strcat('WeightedT_',sbasicCoeffs.labels{2});
        else
            out_dir = sbasicCoeffs.labels{2};
        end

        if tw
            out_dir = strcat(twl,out_dir);

        end


        if isROI
            out_dir = strcat(ROIname,out_dir);
        end

        for npreds= 3:length(sbasicCoeffs.labels)
            out_dir = strcat(out_dir,'_',sbasicCoeffs.labels{npreds});
        end



        out_dir = strcat(out_dir,'_',num2str(thresh));
        mkdir(strcat(PATH, 'Data/Export/', out_dir));

        for subj = 1:length(vps)
            % get the relevant data
            for dataset=1:length(segs)
                clear coeffs
                SEs=[];
                whichdata=segs(dataset);
                if strcmp(whichdata, "r")
                    TIME    = -2000:2:798;

                    BFS     = DATR2;

                    if doViewing
                        looksLeft=rFixationLoc./2;%rlooksLeft;
                    end

                    timeLimit=[-1000 0];
                elseif strcmp(whichdata, "p")
                    TIME    = -2000:2:798;

                    BFS     = DATR2;

                    if doViewing
                        looksLeft=rFixationLoc./2;%rlooksLeft;
                    end
                    timeLimit=[0 500];
                else
                    TIME    = -200:2:3998;

                    BFS     = DAT2;

                    if doViewing
                        looksLeft=FixationLoc./2;%slooksLeft;
                    end
                    timeLimit=[0 1000];
                end


                if isROI

                    if isEOG
                        BFS = nan(1,size(looksLeft,1), size(looksLeft,2));
                        BFS(1,:,:) = abs(looksLeft);
                    else

                        BFS = mean(BFS(ROIelecs, :,:),1);
                    end
                end
                % Matt's script line 1017 - 1090 something
                % we do exclude trials in which participants didn't respond
                fprintf('Processing subject number %d %s-locked\n', vps(subj), whichdata)


                %%%%%%%%%%%%% preparing data for regression %%%%%%%%%%%%%%%%%%

                data = BFS(:,:,allSubDataTable.SubNum==vps(subj) & allSubDataTable.Choice1~=-1);
                bSubData = allSubDataTable(allSubDataTable.SubNum==vps(subj) & allSubDataTable.Choice1~=-1,:);
                if doViewing
                    blooksLeft = looksLeft(:, allSubDataTable.SubNum==vps(subj) & allSubDataTable.Choice1~=-1);
                end

                % Matt actually downsamples
                downSampData=data(:,inc:inc:end,:);
                downSampTimes=TIME(inc:inc:end);
                if doViewing
                    downsampblooksLeft = blooksLeft(inc:inc:end,:);
                end

                if isempty(timeLimit)
                    timeLimit=[min(downSampTimes)-1 max(downSampTimes)+1];
                    subtractor=0;
                else
                    setsubtractor=1;
                end

                timesToLookAt=find(downSampTimes>timeLimit(1) & downSampTimes < timeLimit(2));

                if setsubtractor
                    subtractor=min(timesToLookAt)-1;
                end


                if usefactorscores

                    xMat=meanCentX([ones(size(bSubData,1),1),scoretable{suballSubDataTable.SubNum==vps(subj),1},scoretable{suballSubDataTable.SubNum==vps(subj),2}]);

                else
                    % get subject data, set up regressor matrix, center regressors...

                    if ModsToRun(grabMod) == 1
                        % 1) SUBEVAL MODEL
                        % sbasicCoeffs.labels={'intercept', 'Liking', 'Anxiety', 'Confidence'};
                        xMat=meanCentX([ones(size(bSubData,1),1),  bSubData.Liking/5, bSubData.Anxious/5, bSubData.Confident/5]);

                    elseif  ModsToRun(grabMod) == 2
                        % 2) CHOSEN UNCHOSEN MODEL
                        % sbasicCoeffs.labels={'intercept', 'Chosen', 'Unchosen' }; %
                        xMat=meanCentX([ones(size(bSubData,1),1),  bSubData.ChosenV/10, bSubData.UnchosenV/10]);
                    elseif  ModsToRun(grabMod) == 3
                        % 3) OV RV MODEL
                        % sbasicCoeffs.labels={'intercept', 'VD', 'OV'};
                        xMat=meanCentX([ones(size(bSubData,1),1), bSubData.MaxvMinBid/10, bSubData.AvBid/10]);

                    end

                end
                % loop through channels
                for channel = 1:size(BFS,1)

                    %% loop through time points
                    for t = timesToLookAt
                        % run regression and store coefficients , bSubData.Vlminr/10

                        Y = squeeze(downSampData(channel, t,:));

                        [coeffs(channel,t-subtractor,:), currCIs]  = regress(Y,[xMat]);
                        SEs(channel,t-subtractor,:) =  (currCIs(:,2)- currCIs(:,1))/ 3.92;

                    end
                end
                if strcmp(whichdata, "r")
                    rbasicCoeffs.data(:,:,:,subj)=coeffs;
                    rbasicCoeffs.SEs(:,:,:,subj)=SEs;
                elseif strcmp(whichdata, "p")
                    prbasicCoeffs.data(:,:,:,subj)=coeffs;
                    prbasicCoeffs.SEs(:,:,:,subj)=SEs;
                else
                    sbasicCoeffs.data(:,:,:,subj)=coeffs;
                    sbasicCoeffs.SEs(:,:,:,subj)=SEs;
                end

            end
        end

        save(strcat(PATH, 'Data/Export/', out_dir, '/regcoeffs.mat'), 'sbasicCoeffs', 'rbasicCoeffs', 'prbasicCoeffs')

        %% Cluster-based multiple comparisons corrected hypothesis testing
        % The idea here is to identify "clusters" of activation in the space of electrodes
        % and time.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%         Do significance testing      %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fn=fullfile(strcat(PATH, 'Data/Export'), 'ROI_data.mat');
        regMatsForROI={'sbasicCoeffs','rbasicCoeffs', 'prbasicCoeffs'}; %
        if isROI
            connectionMat=1;
        end
        mkdir(strcat(figDir, out_dir));

        if doSigTesting
            % GET ALL ROIs
            clear allROIs
            for i = 1: length(regMatsForROI)
                eval(['regDat=' regMatsForROI{i} '.data;'])
                eval(['regLabels=' regMatsForROI{i} '.labels;'])
                if weightedT
                    eval(['SEDat=' regMatsForROI{i} '.SEs;'])
                    regDat = regDat.*(1./SEDat);
                end
                whichdata=regMatsForROI{i}(1);
                if strcmp(whichdata, "r")
                    TIME    = -2000:2:798;
                    timeLimit=[-1000 0];
                    %timeLimit=[0 500];
                elseif strcmp(whichdata, "p")
                    TIME    = -2000:2:798;
                    timeLimit=[0 500];

                else
                    TIME    = -200:2:3998;
                    timeLimit=[0 1000];
                end
                downSampTimes=TIME(inc:inc:end);


                if isempty(timeLimit)
                    timeLimit=[min(downSampTimes)-1 max(downSampTimes)+1];
                    subtractor=0;
                else
                    subtractor=min(timesToLookAt)-1;
                end

                timesToLookAt=find(downSampTimes>timeLimit(1) & downSampTimes < timeLimit(2));
                timeLabels=downSampTimes(timesToLookAt);
                for preds = 1:size(regDat, 3) % loop through predictors?



                    EEG_dat=permute(regDat(:,:,preds,:), [4, 1, 2, 3]);
                    % get clusters, cluster sizes, cluster masses for positive clusters (ie
                    % p<.05 in a one tailed positive test):
                    % needs function getEEG_clusterSize
                    % plot and save time course of all predictors


                    posClusterInfo=getEEG_clusterSize_clean(EEG_dat, connectionMat, thresh, 'right');

                    % get clusters, cluster sizes, cluster masses for negative clusters (ie
                    % p<.05 in a one tailed negative test):
                    negClusterInfo=getEEG_clusterSize_clean(EEG_dat, connectionMat, thresh, 'left');



                    % OK... now we have our statistics of interest... we just need a null
                    % distribution to compare them to. Lets run through a loop that:
                    % 1) flips the signs of the data for each subject according to a fair coin
                    % toss.
                    % 2) find the maximum statistic values that you see in the entire dataset
                    % from these "permuted" datasets.
                    % Output:




                    clusterInfo=posClusterInfo;


                    for k=1:numPerm
                        % create a subject length array containing randomly assigned -1's and
                        % 1s:
                        permArray=ones(size(EEG_dat, 1), 1);
                        permArray(logical(binornd(1, .5, size(EEG_dat, 1), 1)))=-1;

                        % multiply each subject timeseries by the -1 or 1 assigned randomly on
                        % this trial
                        sz=size(EEG_dat);
                        permMat=repmat(permArray, [1 sz(2:end)]);

                        % get cluster statistics for permuted dataset:
                        permClusterInfo=getEEG_clusterSize_clean(EEG_dat.*permMat, connectionMat, thresh, 'right');

                        % store the maximum of the statistics, for use in null distribution:
                        maxSize(k)=(max(permClusterInfo.clustSizeMap(:)));
                        maxWt(k)=(max(permClusterInfo.clustWtMap(:)));
                    end


                    % For a two tailed test, find the the minimum cluster statistics necessary
                    % to beat 97.5% of the null distribution
                    % Based on "mass"
                    % I like this statistic better!!!
                    massTh = prctile(maxWt,  97.5 );

                    figure()
                    set(gcf,'Color', [1 1 1])
                    hold on
                    histogram(maxWt)
                    ylimV = get(gca,'ylim');
                    line([massTh massTh],ylimV,'color',[1 0 0], 'linewidth',1)


                    gPos=unique(posClusterInfo.ID_map(posClusterInfo.clustWtMap	>massTh)); % get IDs of significant positive clusters
                    gNeg=unique(negClusterInfo.ID_map(negClusterInfo.clustWtMap	>massTh)); % get IDs of significant negative clusters
                    threshMap=zeros(size(negClusterInfo.tMap));
                    threshMap(posClusterInfo.clustWtMap	>massTh)=posClusterInfo.tMap(posClusterInfo.clustWtMap	>massTh);
                    threshMap(negClusterInfo.clustWtMap	>massTh)=negClusterInfo.tMap(negClusterInfo.clustWtMap	>massTh);
                    nn=[regMatsForROI{i}, '', regLabels{preds}];

                    if ~isROI
                        figure()
                        tint=[1:10];
                        stepsize=10;
                        for tp=1:(size(EEG_dat,3)/stepsize)
                            subplot(5,16,tp)
                            topoplot(mean(EEG_dat(:,:,tint),1),chanlocs,'electrodes','on','plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD);
                            xlim([-0.7 0.7])
                            ylim([-0.7 0.7])
                            title([num2str(timeLabels(tint(1))),':', num2str(timeLabels(tint(end)))]);
                            tint=tint+stepsize;

                        end
                        set(gcf,'units','normalized','position',[.1 .05 .85 .85])
                        [ax,h]=subtitle([nn]);
                        set(h,'Position',[.5 1.05 .5], 'fontsize', 14)

                        fn=[nn];
                        saveas(gcf, fullfile(strcat(figDir, out_dir), fn), 'epsc2')
                        close all

                    end
                    % clusterInfo: Structure containing the following fields:

                    % clustSizeMap: Electrode X Timepoints map of clusters, each labeled with
                    % the number of pairs in that cluster

                    % ID_map: Electrode X Timepoints map of clusters, each labeled with a
                    % unique identifier. Unfortunately, do to my new and improved algorith...
                    % the IDs are not necesarily contigous... so some IDs might be skipped...
                    % but the non-zero IDs will always be unique.

                    % clustWtMap: Electrode X Timepoints map of clusters, each labeled with
                    % the "cluster mass." Cluster mass =  mean|(t-stat)|*pairs in cluster.

                    % tMap: T-statistic from t-test for each electrode*timepoint pair.

                    % STORE ROIs
                    % we'll want:
                    % 1) map --
                    % sign(direction) of effect


                    clear ROIs
                    ROIs.clustThresh = thresh;
                    ROIs.timeLabels = timeLabels;
                    try
                        ROIs.sign=[ones(size(gPos)); -ones(size(gNeg))];
                    catch
                        ROIs.sign=[ones(size(gPos')); -ones(size(gNeg'))]';

                    end
                    k =1;
                    while k <=length(gPos)
                        ROIs.maps(:,:,k)=posClusterInfo.ID_map==gPos(k);
                        ROIs.Tmaps(:,:,k)=ROIs.maps(:,:,k).*abs(posClusterInfo.tMap);
                        % added 03/07/2018 to get size of the cluster
                        ROIs.clustSize(k)= unique(posClusterInfo.clustSizeMap(posClusterInfo.ID_map==gPos(k)));
                        ROIs.wt(k)= unique(posClusterInfo.clustWtMap(posClusterInfo.ID_map==gPos(k)));
                        figure()
                        set(gcf,'Color', [1 1 1])
                        hold on
                        histogram(maxWt)
                        ylimV = get(gca,'ylim');
                        line([massTh massTh],ylimV,'color',[1 0 0], 'linewidth',1)
                        line([ROIs.wt(k) ROIs.wt(k)],ylimV,'color',[1 1 0], 'linewidth',1)
                        % added 07/12/2018 to save out p-values for significant
                        % positive clusters
                        ROIs.p(k) = (1 - (sum(maxWt < ROIs.wt(k)))/length(maxWt)) * 2;
                        % added 03/07/2018 to get involved electrodes and time points
                        clear roielecs
                        clear roitimes
                        [roielecs, roitimes]= find(posClusterInfo.ID_map==gPos(k));
                        ROIs.roielecs(k)={unique(roielecs)};
                        ROIs.roitimes(k)={timeLabels(unique(roitimes))};
                        k=k+1;
                    end

                    while k <= length(gPos) +length(gNeg)
                        ROIs.maps(:,:,k)=negClusterInfo.ID_map==gNeg(k-length(gPos));
                        ROIs.Tmaps(:,:,k)=ROIs.maps(:,:,k).*abs(negClusterInfo.tMap);
                        % added 03/07/2018 to get size of the cluster
                        ROIs.clustSize(k)= unique(negClusterInfo.clustSizeMap(negClusterInfo.ID_map==gNeg(k-length(gPos))));
                        ROIs.wt(k)= unique(negClusterInfo.clustWtMap(negClusterInfo.ID_map==gNeg(k-length(gPos))));
                        % added 07/12/2018 to save out p-values for significant
                        % negative clusters
                        ROIs.p(k) = (1 - (sum(maxWt < ROIs.wt(k)))/length(maxWt)) * 2;
                        % added 03/07/2018 to get involved electrodes and time points
                        figure()
                        set(gcf,'Color', [1 1 1])
                        hold on
                        histogram(maxWt)
                        ylimV = get(gca,'ylim');
                        line([massTh massTh],ylimV,'color',[1 0 0], 'linewidth',1)
                        line([ROIs.wt(k) ROIs.wt(k)],ylimV,'color',[1 1 0], 'linewidth',1)
                        clear roielecs
                        clear roitimes
                        [roielecs, roitimes]= find(negClusterInfo.ID_map==gNeg(k-length(gPos)));
                        ROIs.roielecs(k)={unique(roielecs)};
                        ROIs.roitimes(k)={timeLabels(unique(roitimes))};
                        k=k+1;
                    end

                    for k = 1:length(ROIs.sign)

                        clustTs=ROIs.maps(:,:,k).*abs(posClusterInfo.tMap); % this is only getting posCluster info...
                        % apparently all t-values are in both maps (they're
                        % identical, so this is fine.
                        [ROIs.peakChannel(k),J] = find(clustTs==max(clustTs(:)));

                        ROIs.peakTime(k)=timeLabels(J);
                    end


                    if ~ isROI
                        % MAKE basic plots:
                        figure();
                        subplot(2, 1,2)

                        cLim=[8];
                        hold on
                        imagesc(timeLabels, 1:length(chanlocs), posClusterInfo.tMap)
                        set(gca, 'clim', [-cLim, cLim])
                        colorbar
                        xlim([min(timeLabels), max(timeLabels)])
                        ylim([0, nchans])
                        set(gca, 'box', 'off', 'ytickLabels', '')
                    end

                    if isfield(ROIs, 'peakTime')
                        plot(ROIs.peakTime,  ROIs.peakChannel, 'ok', 'markerSize', 6, 'markerFaceColor', 'none', 'markerEdgeColor', 'w', 'lineWidth', 1)
                    end

                    %             eventTimes=[0, 500, 1500, 2000];
                    %             for event=1:length(eventTimes)
                    %                 plot([eventTimes(event), eventTimes(event)], [1, 64], '-k')
                    %             end


                    fn=[nn, 'just_time.eps'];
                    saveas(gcf, fullfile(strcat(figDir, out_dir),fn), 'epsc2')
                    close all


                    if ~isROI
                        if isfield(ROIs, 'peakTime')
                            for k = 1:length(ROIs.peakTime)
                                figure
                                topoplot(posClusterInfo.tMap(:,timeLabels==ROIs.peakTime(k)),chanlocs,'electrodes','off','plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {ROIs.roielecs{k},'*','k', 3,1});
                                xlim([-0.7 0.7])
                                ylim([-0.7 0.7])
                                title([nn, ' ' num2str(ROIs.peakTime(k))]);
                                fn=[nn, ' ' num2str(ROIs.peakTime(k))];
                                saveas(gcf, fullfile(strcat(figDir, out_dir), fn), 'epsc2')
                                close all
                            end

                        end
                    end

                    % MAKE plot of individual subject ROI values:


                    % get subject ROI coefficients:

                    if isfield(ROIs, 'maps')



                        subCoefMat=nan(size(EEG_dat, 1), size(ROIs.maps,3));
                        for k = 1:size(EEG_dat, 1);
                            subDat=squeeze(EEG_dat(k,:,:));
                            for z = 1:size(ROIs.maps,3)
                                subCoefMat(k,z)=nanmean(subDat(ROIs.maps(:,:,z)));
                            end
                        end

                        ROIs.subCoeffs=subCoefMat;

                        Scale=max(abs(subCoefMat(:)));
                        xJit = smartJitter(subCoefMat,.1,.1);

                        ll=size(subCoefMat, 1);
                        plot([0, size(subCoefMat, 2)], [0 0], '--k')

                        nPlot=size(subCoefMat, 2);

                        hold on
                        plot([0, nPlot+.5], [0, 0], '--k')

                        for k = 1:nPlot
                            plot(ones(ll, 1).*k+xJit(:,k) , subCoefMat(:,k), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(k,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
                        end
                        ylim([-Scale, Scale]);
                        xlim([.5, size(subCoefMat, 2)+.5]);
                        ylabel('Coefficient')
                        set(gca, 'xticklabel', {}, 'box', 'off')

                        fn=fullfile(strcat(figDir, out_dir), [nn, '_indSubCoeffs.eps']);
                        saveas(gcf, fn, 'epsc2')
                        close all
                    end
                    eval(['allROIs.', nn, '=ROIs;']);
                end
            end
            fn=fullfile(strcat(PATH, 'Data/Export/', out_dir), 'ROI_data.mat');

            save(fn, 'allROIs');
        else
            load(fn);
        end

    end
end




%% here's how I generated DAT2 and DATR2
% Set everything prior to stimulus onset to nan for S-locked trials
RTms = allSubDataTable.RTeval1*1000;
RTtp = round(RTms/2)+100;
RTtp(isnan(RTtp))=1; % set nan values to one to start nan-ing out from beginning for no-response trials
DAT2 = DAT;

for tidx=1:size(DAT,3)
    DAT2(:,RTtp(tidx):size(DAT2,2),tidx) = nan;

end

% Set everything prior to stimulus onset to nan for R-locked trials
RTms = allSubDataTable.RTeval1*1000;

% + 800 is basically our new reference (it's where the data ends and we go
% backwards from that
% so we add 400 sampling points to our half RT in ms to get the amount of
% sampling points covered and safe. (if RT + post Rt is longer than our
% segment length, we don't have no trial data. Otherwise we need to
% eliminate all info happening in the gap between those.
RTtp = round(RTms/2)+400; % now we have 800, so we need to add 400 sampling points instead
RTtp(isnan(RTtp))=size(DATR,3); % set nan values to one
DATR2 = DATR;

for tidx=1:size(DATR,3) % for all trials
    if RTtp(tidx)<size(DATR2,2) % check whether we're safe or need to do something
        DATR2(:,1:(size(DATR2,2)-RTtp(tidx)),tidx) = nan; % if we need to, nan
        % out everything that happens before the trial onset (1st cell to where the trial starts:
        % end of segment - RT + post resp.
    end
end

% before excluding EOG channels, save out the data that had them



oDATR2 = DATR2; % save original data
oDAT2 = DAT2;
DATR2 = DATR2([1:60,65],:,:); % get rid of EOG chans
DAT2 = DAT2([1:60,65],:,:);
save(sprintf('%sData/Export/DAT2new_wEOG.mat', PATH), 'oDAT2', '-v7.3');
save(sprintf('%sData/Export/DATR2new_wEOG.mat', PATH), 'oDATR2', '-v7.3');
save(sprintf('%sData/Export/DAT2new.mat', PATH), 'DAT2', '-v7.3');
save(sprintf('%sData/Export/DATR2new.mat', PATH), 'DATR2', '-v7.3');

%% Generate factor scores
% 1b) exclude trials with no response!5028
% get factor scores

% indicators2=[allSubDataTable.ChosenV(allSubDataTable.Choice1~=-1 & allSubDataTable.SubNum~=5028)/10,allSubDataTable.UnchosenV(allSubDataTable.Choice1~=-1 & allSubDataTable.SubNum~=5028)/10, allSubDataTable.MaxvMinBid(allSubDataTable.Choice1~=-1 & allSubDataTable.SubNum~=5028)/10, allSubDataTable.AvBid(allSubDataTable.Choice1~=-1 & allSubDataTable.SubNum~=5028)/10, allSubDataTable.AvBidSetSal(allSubDataTable.Choice1~=-1 & allSubDataTable.SubNum~=5028)/5,allSubDataTable.Anxious(allSubDataTable.Choice1~=-1 & allSubDataTable.SubNum~=5028)/5, allSubDataTable.Liking(allSubDataTable.Choice1~=-1 & allSubDataTable.SubNum~=5028)/5, allSubDataTable.Confident(allSubDataTable.Choice1~=-1 & allSubDataTable.SubNum~=5028)/5];
% [coeff,score,latent,tsquared,explained,mu] = pca(indicators2,'NumComponents',2); % change back to 3
% [L1,T] = rotatefactors(coeff);
% score= indicators2*L1;
% scoretable = array2table(score);
% scoretable.Properties.VariableNames = {'Valuation','Choice'};
% save(sprintf('%sData/Export/PCAfactorScoreTable.mat', PATH), 'scoretable')

%% GET CONNECTION MATRIX specifying which electrodes are connected to which:
% % % R: why 40? what's the measure? mm
% threshold distance for "connected" electrodes
recreateConnectionMat=0;
if  recreateConnectionMat
    connectThresh=40;
    % if we've already got one made, load it.
    % otherwise, create one from scratch:
    clear relDist
    % Then get the XYZ coordinates for each channel:
    allLocas=[[chanlocs.X]; [chanlocs.Y]; [chanlocs.Z]]' ;
    % Loop through the channels and get the distance between that channel and
    % all other channels.
    for i = 1:length(allLocas)
        relDist(:,i)=sqrt(sum((allLocas-repmat(allLocas(i,:), length(allLocas), 1)).^2, 2));
    end

    % Set a threshold on distance... and mark channels that fall within that
    % threshold:
    connectionMat=relDist<connectThresh;
    connectionMat=connectionMat-eye(length(connectionMat)); % eye makes identity matrix

    % CHECK OUT CONNECTIONS:
    %imagesc(connectionMat);
    %close all

    % check average number of connections:
    %numCons=nanmean(sum(connectionMat));

    % Look at connectivity weights across scalp, alla AC:
    % figure;
    % for e = 1:64
    % subplot(8,8,e)
    % topoplot(connectionMat(e,:),EEG.chanlocs)
    % end
    % saveas(gcf, 'connectionMatFig.eps', 'epsc2')
    % close all
    save( strcat(PATH, 'Data/Export/connectionMat.mat'), 'connectionMat')
end

%% compute viewing direction (based on original data)
% I decided to clean up the EOG data a bit because there are some quite
% weird values in there...
% doViewing=1;
% if doViewing
% %     load(sprintf('%sData/Export/DAT2new_wEOG.mat', PATH))
% %     load(sprintf('%sData/Export/DATR2new_wEOG.mat', PATH))
%     for n= 1: length(vpss)
%         vpidx = find(allSubDataTable.SubNum==vpss(n));
%         %MAXAMPD(vpidx) = nanmedian(max(abs(oDAT2(LO1,:,vpidx)-oDAT2(LO2,:,vpidx))));
%         MAXDIST=reshape(max(abs(oDAT2(LO1,:,vpidx)-oDAT2(LO2,:,vpidx))), [1 length(vpidx)]);
%         cutp=prctile(MAXDIST,95);
%         medV = nanmedian(MAXDIST);
%         sdcut= 3*std(MAXDIST, 'omitnan');
%         [likr,centersr]= ksdensity(MAXDIST, 'Support', 'positive');
%         figure()
%         hold on
%         histfit(MAXDIST, 50, 'kernel')
%         plot([cutp cutp], ylim, '--r')
%         plot([medV medV], ylim, '--b')
%         plot([medV+sdcut medV+sdcut], ylim, '--b')
%
%         isoutlier = abs(oDAT2(LO1,:,vpidx)-oDAT2(LO2,:,vpidx)) > medV+sdcut;
%         isoutlier(2,:,:) = isoutlier;
%         %ksdensity(MAXDIST, 'Support', 'positive')
%         tmpdat = oDAT2([LO1 LO2],:,vpidx);
%         tmpdat(isoutlier) = nan;
%         oDAT2([LO1 LO2],:,vpidx)=tmpdat;
%
%         MAXDIST=reshape(max(abs(oDAT2(LO1,:,vpidx)-oDAT2(LO2,:,vpidx))), [1 length(vpidx)]);
%         cutp=prctile(MAXDIST,95);
%         medV = nanmedian(MAXDIST);
%         sdcut= 3*std(MAXDIST, 'omitnan');
%         [likr,centersr]= ksdensity(MAXDIST, 'Support', 'positive');
%         figure()
%         hold on
%         title('check')
%         histfit(MAXDIST, 50, 'kernel')
%         plot([cutp cutp], ylim, '--r')
%         plot([medV medV], ylim, '--b')
%         plot([medV+sdcut medV+sdcut], ylim, '--b')
%
%         %     MINAMPD(vpidx)= nanmedian(min(abs(oDAT2(LO1,:,vpidx)-oDAT2(LO2,:,vpidx))));
%         %     MED(vpidx)=nanmedian(abs(oDAT2(LO1,:,vpidx)-oDAT2(LO2,:,vpidx)));
%     end
%     save(sprintf('%sData/Export/DAT2new_wEOG_tidied_up.mat', PATH), 'oDAT2', '-v7.3');
%     save(sprintf('%sData/Export/DATR2new_wEOG_tidied_up.mat', PATH), 'oDATR2', '-v7.3');
%     close all
%
%     slooksLeftnn = reshape((oDAT2(LO1,:,:)-oDAT2(LO2,:,:)), [size(oDAT2,2) size(oDAT2,3)]);
%
%     % slocked
%     slooksLeft = reshape((oDAT2(LO1,:,:)-oDAT2(LO2,:,:))./(2*(abs(oDAT2(LO1,:,:))+ abs(oDAT2(LO2,:,:)))), [size(oDAT2,2) size(oDAT2,3)]);
%
%     % r-locked
%     rlooksLeft = reshape((oDATR2(LO1,:,:)-oDATR2(LO2,:,:))./(2*(abs(oDATR2(LO1,:,:))+ abs(oDATR2(LO2,:,:)))), [size(oDATR2,2) size(oDATR2,3)]);
%
%
% % slooksLeft without saccades
%
%     % first find stable intervals:
%
%     isfixation = zeros(size(slooksLeft));
%     for n=1:size(slooksLeft,1)-1
%         isfixation(n,:)=(slooksLeft(n+1,:)==slooksLeft(n,:)& ~isnan(slooksLeft(n,:)));
%
%     end
%
%
%
% % set saccades to zero
%
%     FixationLoc = slooksLeft*2;
%     FixationLoc(~isfixation)=0;
%     isfixation(isnan(slooksLeft))= nan;
%
%
% % same thing for response locked
%
%     risfixation = zeros(size(rlooksLeft));
%     for n=1:size(rlooksLeft,1)-1
%         risfixation(n,:)=(rlooksLeft(n+1,:)==rlooksLeft(n,:)& ~isnan(rlooksLeft(n,:)));
%
%     end
%
% % set saccades to zero for resp-locked
%     rFixationLoc = rlooksLeft*2;
%     rFixationLoc(~risfixation)=0;
%     risfixation(isnan(rlooksLeft))= nan;
%
%
%     save(sprintf('%sData/Export/RnSFixationLoc.mat', PATH), 'rFixationLoc', 'FixationLoc', '-v7.3');
%
%     % check
%     CAmps = slooksLeftnn;
%     CAmps(FixationLoc==0)= 0;
%
% % naning out non-trial
%     RTms = allSubDataTable.RTeval1*1000;
%     RTtp = round(RTms/2)+100;
%     RTtp(isnan(RTtp))=1; % set nan values to one
%     FixLocALL = FixationLoc;
%
%     for tidx=1:size(FixationLoc,2)
%         FixationLoc(RTtp(tidx):size(FixationLoc,1),tidx) = nan;
%
%     end
% % retrieve oscillations per choice
%     isoscillation = zeros(size(FixationLoc));
%     lastFix = nan(1,size(FixationLoc,2));
%     lastFix(FixationLoc(1,:)~=0 & ~isnan(FixationLoc(1,:))) = FixationLoc(1, FixationLoc(1,:)~=0 & ~isnan(FixationLoc(1,:))); % get fixations at first time point
%     tmpFixationLoc = FixationLoc;
%     for n=2:size(FixationLoc,1)
%
%         tmpFixationLoc(n, FixationLoc(n,:)==0) = lastFix( FixationLoc(n,:)==0);
%         isoscillation(n,:)=(tmpFixationLoc(n,:)~=tmpFixationLoc(n-1,:) & FixationLoc(n,:)~=0 & ~isnan(FixationLoc(n,:))); % get all the values where stuff actually differs from one another
%         % now overwrite the values that are not actually sign switches
%
%         lastFix(FixationLoc(n,:)~=0 & ~isnan(FixationLoc(n,:))) = FixationLoc(n, FixationLoc(n,:)~=0 & ~isnan(FixationLoc(n,:)));
%
%     end
%
%     noscillations=sum(isoscillation(101:end,:),1);
%
%     allSubDataTable.nOsc= noscillations';
% %
% %     figure()
% %     hold on
% %     %title('check')
% %     histfit(noscillations)
% % add fixation info to table for export
%     allSubDataTable.fixLeft= sum(FixationLoc(101:end,:)==1)'*2;
%     allSubDataTable.fixRight= sum(FixationLoc(101:end,:)== -1)'*2;
%     allSubDataTable.relRight= (sum(FixationLoc(101:end,:)== -1)' - sum(FixationLoc(101:end,:)==1)')./(sum(FixationLoc(101:end,:)== -1)' + sum(FixationLoc(101:end,:)==1)');
%
% % get fixation-location within the first 100 and 300 ms
%     allSubDataTable.firstfix= mean(FixationLoc(192:369,:), 1)';
%     allSubDataTable.firstfixbin = ones(size(allSubDataTable.firstfix));
%     allSubDataTable.firstfixbin(allSubDataTable.firstfix > 0) = -1;
% % save stuff out
%     EOGallSubDataTable = allSubDataTable;
%     writetable(EOGallSubDataTable, sprintf('%sData/Export/EOGallSubDataTable.xls', PATH));
%
% % export FixLoc
%     FixLocExp0 = FixationLoc(1:800,:)';
%
%     FixLocExp = nan(size(FixLocExp0,1), 150);
%     for s= 1:150
%         tw = 1+((s-1)*5):(s*5);
%         FixLocExp(:,s) = nanmean(FixLocExp0(:,tw),2);
%     end
%
%     FixLocExp = array2table(FixLocExp);
%
%     clear('prTFh')
%     for s= 1:150
%         prTFh{s}= strcat('tw',num2str(abs(-200+(s-1)*10)), num2str(abs(-200+(s)*10)));
%
%     end
%
%     FixLocExp.Properties.VariableNames = prTFh;
%     writetable(FixLocExp, sprintf('%s%sFixTime.xls', PATH, EXPORTPATH));
%
% % do just unsigned fixation
% % export FixLoc
%     FixLocExp0 = abs(FixationLoc(1:800,:)');
%
%     FixLocExp = nan(size(FixLocExp0,1), 150);
%     for s= 1:150
%         tw = 1+((s-1)*5):(s*5);
%         FixLocExp(:,s) = nanmean(FixLocExp0(:,tw),2);
%     end
%
%     FixLocExp = array2table(FixLocExp);
%
%     clear('prTFh')
%     for s= 1:150
%         prTFh{s}= strcat('tw',num2str(abs(-200+(s-1)*10)), num2str(abs(-200+(s)*10)));
%
%     end
%
%     FixLocExp.Properties.VariableNames = prTFh;
%
%     writetable(FixLocExp, sprintf('%s%sFixTimeabs.xls', PATH, EXPORTPATH));
% % export r-locked fixation
%
% % export FixLoc
%     FixLocExp0 = rFixationLoc(end-900:end-101,:)';
%
%     FixLocExp = nan(size(FixLocExp0,1), 150);
%     for s= 1:150
%         tw = 1+((s-1)*5):(s*5);
%         FixLocExp(:,s) = nanmean(FixLocExp0(:,tw),2);
%     end
%
%     FixLocExp = array2table(FixLocExp);
%
% %
%     clear('prTFh')
%     for s= 1:150
%         prTFh{s}= strcat('tw',num2str(abs(-1500+(s-1)*10)), num2str(abs(-1500+(s)*10)));
%
%     end
%
%     FixLocExp.Properties.VariableNames = prTFh;
%
%     writetable(FixLocExp, sprintf('%s%sFixTimer.xls', PATH, EXPORTPATH));
%
% % same: absolute fixation
% % export FixLoc
%     FixLocExp0 = abs(rFixationLoc(end-900:end-101,:)');
%
%     FixLocExp = nan(size(FixLocExp0,1), 150);
%     for s= 1:150
%         tw = 1+((s-1)*5):(s*5);
%         FixLocExp(:,s) = nanmean(FixLocExp0(:,tw),2);
%     end
%
%     FixLocExp = array2table(FixLocExp);
%
%     clear('prTFh')
%     for s= 1:150
%         prTFh{s}= strcat('tw',num2str(abs(-1500+(s-1)*10)), num2str(abs(-1500+(s)*10)));
%
%     end
%
%     FixLocExp.Properties.VariableNames = prTFh;
%
%     writetable(FixLocExp, sprintf('%s%sFixTimerabs.xls', PATH, EXPORTPATH));
%
% end

%%
function [ax,h]=subtitle(text)
%
%Centers a title over a group of subplots.
%Returns a handle to the title and the handle to an axis.
% [ax,h]=subtitle(text)
%           returns handles to both the axis and the title.
% ax=subtitle(text)
%           returns a handle to the axis only.
ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on')
title(text);
if (nargout < 2)
    return
end
h=get(ax,'Title');

end



