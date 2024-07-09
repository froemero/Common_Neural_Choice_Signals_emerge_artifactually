%% To do mass-univariate ala Matt (Nassar)
% script adapted from Matt's cannon task script
% if anything breaks, contact romy.froemer@googlemail.com

% load data
PATH = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Analyses/Matlab/UCAP/';
figDir= strcat(PATH, 'Figures/');

%% beta data see below for how this is generated
concatufresults=0; % set to 1 if you need to generate the relevant matrix and run that below

load(strcat(PATH, 'Second_level/Export/regcoeffs.mat')) % this loads stimulus and response-locked concatenated uf results plus chanlocs
% load this code to make it is at the end of the script
load(strcat(PATH, 'Second_level/connectionMat.mat'))
%noEOG=1; % with or without EOG channels

%% get channel locations and channel numbers (for plotting)
load(strcat(PATH, 'Second_level/Export/chanlocs.mat'))
nchans = size(chanlocs,2);
for n = 1:nchans
    
    expression = [chanlocs(n).labels '=' sprintf('%d',n) ';'];
    eval(expression)
end
% %% in this script we don't use EOG, so we exclude those channels
% nchans=61; % reset number of channels
% chanlocs= chanlocs([1:60,65]); % adjust chanlocs
% Cz=61; % adjust channumber for Cz



%% generate path to helper functions
addpath(strcat('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/', 'helpfunctions'))

% contains:
% meanCentx - for centering predictors
% nanzscore - the name says it
% zScoreX   - I assume that's being used by some of the other functions?
% getEEG_clusterSize - to obtain positive and negative clusters during
% permutation test

%%
% set parameters

whichdata = "s"; % "r" or "s"
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
doSigTesting=1;
isROI=0;
inc = 1; % we don't downsample
weightedT=0; % can't do this here, because I don't have the SEs
numPerm=1000;
ocbColors=[0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167]./256;
cbColors= [ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors; ocbColors];

%% Here we loop through a vector that spacifies which analyses should be run

AllThreshs= [.001,.005];%

for currthresh= 1:length(AllThreshs) % loop through thresholds I want to test
    thresh=AllThreshs(currthresh);
    isEOG =0;
    
    %% Cluster-based multiple comparisons corrected hypothesis testing
    % The idea here is to identify "clusters" of activation in the space of electrodes
    % and time.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         Do significance testing      %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fn=fullfile(strcat(PATH, 'Second_level/Export'), 'ROI_data.mat');
    regMatsForROI={'sbasicCoeffs','rbasicCoeffs', 'prbasicCoeffs'}; %
    
    out_dir = sbasicCoeffs.labels{2};
    for npreds= 3:length(sbasicCoeffs.labels)
        out_dir = strcat(out_dir,'_',sbasicCoeffs.labels{npreds});
    end
    out_dir = strcat(out_dir,'_',num2str(thresh));
    
    mkdir(strcat(PATH, 'Second_level/Export/', out_dir));
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
                TIME    = rbasicCoeffs.times *1000;
                timeLimit=[-1000 200];
            elseif strcmp(whichdata, "p")
                TIME    = prbasicCoeffs.times *1000;
                timeLimit=[-200 998];
            else
                TIME    = sbasicCoeffs.times *1000;
                timeLimit=[-200 998];
            end
            downSampTimes=TIME(inc:inc:end);
            
            
            if isempty(timeLimit)
                timeLimit=[min(downSampTimes) max(downSampTimes)];%timeLimit=[min(downSampTimes)-1 max(downSampTimes)+1];
                subtractor=0;
                
            end
            timesToLookAt=find(downSampTimes>=timeLimit(1) & downSampTimes <= timeLimit(2));
            %             if ~isempty(timeLimit)
            %                 subtractor=min(timesToLookAt)-1;
            %             end
            
            
            timeLabels=downSampTimes(timesToLookAt);
            for preds = 1:size(regDat, 3) % loop through predictors?
                
                
                %% see whether this actually works
                EEG_dat=permute(regDat(:,timesToLookAt,preds,:), [4, 1, 2, 3]);
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
                        try
                            title([num2str(timeLabels(tint(1))),':', num2str(timeLabels(tint(end)))]);
                            tint=tint+stepsize;
                        catch
                            title([num2str(timeLabels(tint(1))),':', num2str(timeLabels(end))]);
                            tint=tint+stepsize;
                        end
                        
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
        fn=fullfile(strcat(PATH, 'Second_level/Export/', out_dir), 'ROI_data.mat');
        
        save(fn, 'allROIs');
    else
        load(fn);
    end
    
end

%% I need to concatenate all subjects' results before I can do the permutation test

if concatufresults
    %% loop through subjects, generate and save
    SOURCEFILES = dir(strcat(PATH, 'UCAP_UF_Results/*.mat')); %all cleaned files
    SUBJECTS = 1:numel(SOURCEFILES);
    clear 'sbasicCoeffs'
    clear 'rbasicCoeffs'
    clear 'prbasicCoeffs'
    
    sbasicCoeffs.labels={'intercept', 'DeviantP', 'Quality', 'DPbyQuality'};
    rbasicCoeffs.labels={'intercept', 'DeviantP', 'Quality', 'DPbyQuality'};
    prbasicCoeffs.labels={'intercept', 'DeviantP', 'Quality', 'DPbyQuality'};
    
    
    
    for subj = 1:length(SUBJECTS)
        
        % load dataset
        load(sprintf('%sUCAP_UF_Results/%s',PATH, SOURCEFILES(subj).name))
        if subj == 1
            
            sbasicCoeffs.times= ufresult.times([find(ufresult.times== -0.2) : find(ufresult.times== 0.998)]);
            rbasicCoeffs.times= ufresult.times([find(ufresult.times== -1) : find(ufresult.times== 0.2)]);
            prbasicCoeffs.times=ufresult.times([find(ufresult.times== -0.2) : find(ufresult.times== 0.998)]);
        end
        % concatenate datasets and save them as separate structures for s, r,
        % and post r
        sbasicCoeffs.data(:,:,:,subj)=ufresult.beta(:,[find(ufresult.times== -0.2) : find(ufresult.times== 0.998)], 1:4);
        
        rbasicCoeffs.data(:,:,:,subj)=ufresult.beta(:,[find(ufresult.times== -1) : find(ufresult.times== 0.2)], 5:8);
        
        prbasicCoeffs.data(:,:,:,subj)=ufresult.beta(:,[find(ufresult.times== -0.2) : find(ufresult.times== 0.998)], 5:8);
        
    end
    chanlocs = ufresult.chanlocs;
    % save
    save(strcat(PATH, 'Second_level/Export/', 'regcoeffs.mat'), 'sbasicCoeffs', 'rbasicCoeffs', 'prbasicCoeffs', 'chanlocs')
    
    clear ufresult
end
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
    save( strcat(PATH, 'Second_level/connectionMat.mat'), 'connectionMat')
end


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



