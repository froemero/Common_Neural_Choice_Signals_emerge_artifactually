
%% Goals:
% generate Figures of MASS vs Unfold results for the 4 datasets we
% reanalyzed (Pisauro, UCAP, Steinemann & Boldt)

% For each we will load in the y-hats to plot curves, the RTs to plot
% histograms, the betas at the relevant electrode (Pz if possible,
% otherwise CPz), to do t-tests across the relevant time-range and see
% whether there's a significant difference (we'll do that for each, unfold and Mass results).

basepath = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Data/';
figpath = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Docu/Reanalysis_Figures/';

addpath('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/helpfunctions');
figsize= [461   119   304   460];
numPerm=1000;

%% Pisauro

FilePath = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Data/Pisauro/';

allFiles = dir('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Data/Pisauro/sub-*_EEG_EEG_events_sub*.mat');

alldata = {};

for subi = 1: size(allFiles,1)

subID = allFiles(subi).name(5:6);

load(sprintf('%s/sub-%s_EEG_EEG_events_sub%s.mat', FilePath, subID, subID))

subRT = [events{1}.rt'; events{2}.rt'];

subdata = table(repmat(subID, length(subRT),1), subRT);

alldata = [alldata; subdata];

clear events
end
%%

%plot RTs

sRTms = round(alldata.subRT*1000);
sRTms = alldata.subRT;
%% cf median split on all RTs
f5=find(sRTms<= quantile(sRTms, 0.5));
f1=find(sRTms>quantile(sRTms, 0.5));



Pisauro = readtable(sprintf('%sPisauro/post_unfold/2023-10-11_bothRT', basepath));

% separate out: stimulus vs response; deconv vs epoched

Pisauro_s_deconv = Pisauro(contains(Pisauro.type,'deconv') & contains(Pisauro.basisname, 'stimulus') & Pisauro.time> -0.2 & Pisauro.time< 1,:);
Pisauro_r_deconv = Pisauro(contains(Pisauro.type,'deconv') & contains(Pisauro.basisname, 'response') & Pisauro.time> -1 & Pisauro.time< 0.2,:);

Pisauro_s_MASS = Pisauro(contains(Pisauro.type,'epoched') & contains(Pisauro.basisname, 'stimulus') & Pisauro.time> -0.2 & Pisauro.time< 1,:);
Pisauro_r_MASS = Pisauro(contains(Pisauro.type,'epoched') & contains(Pisauro.basisname, 'response') & Pisauro.time> -1 & Pisauro.time< 0.2,:);

Pisauro_r_deconv_beta = Pisauro.yhat(contains(Pisauro.rt_split, 'slower=') & contains(Pisauro.type,'deconv') & contains(Pisauro.basisname, 'response') & Pisauro.time> -1 & Pisauro.time< 0.2)...
    - Pisauro.yhat(contains(Pisauro.rt_split, 'faster') & contains(Pisauro.type,'deconv') & contains(Pisauro.basisname, 'response') & Pisauro.time> -1 & Pisauro.time<0.2);

Pisauro_r_MASS_beta = Pisauro.yhat(contains(Pisauro.rt_split, 'slower=') & contains(Pisauro.type,'epoched') & contains(Pisauro.basisname, 'response') & Pisauro.time> -1 & Pisauro.time< 0.2)...
    - Pisauro.yhat(contains(Pisauro.rt_split, 'faster') & contains(Pisauro.type,'epoched') & contains(Pisauro.basisname, 'response') & Pisauro.time> -1 & Pisauro.time< 0.2);

Pisauro_s_deconv_beta = Pisauro.yhat(contains(Pisauro.rt_split, 'slower=') & contains(Pisauro.type,'deconv') & contains(Pisauro.basisname, 'stimulus') & Pisauro.time> -0.2 & Pisauro.time< 1)...
    - Pisauro.yhat(contains(Pisauro.rt_split, 'faster') & contains(Pisauro.type,'deconv') & contains(Pisauro.basisname, 'stimulus') & Pisauro.time> -0.2 & Pisauro.time< 1);

Pisauro_s_MASS_beta = Pisauro.yhat(contains(Pisauro.rt_split, 'slower=') & contains(Pisauro.type,'epoched') & contains(Pisauro.basisname, 'stimulus') & Pisauro.time> -0.2 & Pisauro.time< 1)...
    - Pisauro.yhat(contains(Pisauro.rt_split, 'faster') & contains(Pisauro.type,'epoched') & contains(Pisauro.basisname, 'stimulus') & Pisauro.time> -0.2 & Pisauro.time< 1);


Pisauro_r_time = Pisauro.time(contains(Pisauro.rt_split, 'slower=') & contains(Pisauro.type,'deconv') & contains(Pisauro.basisname, 'response') & Pisauro.time> -1 & Pisauro.time< 0.2);

Pisauro_r_subs = Pisauro.subject(contains(Pisauro.rt_split, 'slower=') & contains(Pisauro.type,'deconv') & contains(Pisauro.basisname, 'response') & Pisauro.time> -1 & Pisauro.time< 0.2);%% goal: compute t-test at all time-points

%%
% sub-goal: make matrix where columns are Subjects and timepoints are
% columns

All_Betas_Pisauro_MASS = nan(length(unique(Pisauro_r_subs)),length(unique(Pisauro_r_time)));
All_Betas_Pisauro_deconv = nan(length(unique(Pisauro_r_subs)),length(unique(Pisauro_r_time)));

All_y_Pisauro_MASS = nan(length(unique(Pisauro_r_subs)),length(unique(Pisauro_r_time)), 2);
All_y_Pisauro_deconv = nan(length(unique(Pisauro_r_subs)),length(unique(Pisauro_r_time)), 2);

tps = unique(Pisauro_r_time);
for tp = 1:length(tps)
    All_Betas_Pisauro_MASS(:,tp) = Pisauro_r_MASS_beta(Pisauro_r_time==tps(tp));
    All_Betas_Pisauro_deconv(:,tp) = Pisauro_r_deconv_beta(Pisauro_r_time==tps(tp));
    All_y_Pisauro_MASS(:,tp, 1) = Pisauro_r_MASS.yhat(contains(Pisauro_r_MASS.rt_split, 'slower=')& Pisauro_r_MASS.time==tps(tp));
    All_y_Pisauro_MASS(:,tp, 2) = Pisauro_r_MASS.yhat(contains(Pisauro_r_MASS.rt_split, 'faster')& Pisauro_r_MASS.time==tps(tp));
    All_y_Pisauro_deconv(:,tp, 1) = Pisauro_r_deconv.yhat(contains(Pisauro_r_deconv.rt_split, 'slower=')& Pisauro_r_deconv.time==tps(tp));
    All_y_Pisauro_deconv(:,tp, 2) = Pisauro_r_deconv.yhat(contains(Pisauro_r_deconv.rt_split, 'faster')& Pisauro_r_deconv.time==tps(tp));
end

%% now do t-test
thresh = 0.05;
[sig_Pisauro_MASS, pM, ~, dM]=ttest(All_Betas_Pisauro_MASS, [], thresh, 'both');
All_t_Pisauro_MASS = dM.tstat;
[sig_Pisauro_deconv, pD, ~, dD]=ttest(All_Betas_Pisauro_deconv, [], thresh, 'both');
All_t_Pisauro_deconv = dD.tstat;

sig_Pisauro_MASS(sig_Pisauro_MASS == 0) = nan;
sig_Pisauro_deconv(sig_Pisauro_deconv == 0) = nan;

sigClusts_Pisauro= [tps(sig_Pisauro_deconv==1), pD(sig_Pisauro_deconv==1)'];

%% permutation correction
timeLabels = tps;
PermStats_Pisauro_MASS = clustercorrection(All_Betas_Pisauro_MASS,[1], thresh, numPerm,timeLabels);
sig_Pisauro_MASS = sum(PermStats_Pisauro_MASS.maps,3);
sig_Pisauro_MASS(sig_Pisauro_MASS == 0) = nan;
PermStats_Pisauro_Deconv = clustercorrection(All_Betas_Pisauro_deconv,[1], thresh, numPerm,timeLabels);
try
sig_Pisauro_deconv = sum(PermStats_Pisauro_Deconv.maps,3);
sig_Pisauro_deconv(sig_Pisauro_deconv == 0) = nan;
catch
sig_Pisauro_deconv = nan(size(sig_Pisauro_MASS));
end



%%
figure()
set(gcf,'Color', [1 1 1])
hold on
Leg1 = plot(tps, All_t_Pisauro_MASS, 'k-', LineWidth=2);
plot(tps, sig_Pisauro_MASS*-8, 'k-', LineWidth=4)
Leg2 = plot(tps, All_t_Pisauro_deconv, '-','color', [0.6 0.6, 0.6], LineWidth=2);
plot(tps, sig_Pisauro_deconv*-7.5, '-','color', [0.6 0.6, 0.6], LineWidth=4);
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'Mass-Univariate', 'Unfold-Deconvolution');% , Leg3, Leg4
set(gca,'TickDir','out', 'Box', 'off')
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14, 'position',[.25,.25,lp(3:4)])%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')

%% let's just plot the conditions:

figure()
set(gcf,'Color', [1 1 1])
subplot(7,1,[1:3])
hold on
plot(tps,mean(All_y_Pisauro_MASS(:,:,1),1) , '-','color', [107/255 0, 0], LineWidth=2)
plot(tps, mean(All_y_Pisauro_MASS(:,:,2),1), '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, sig_Pisauro_MASS*-1.5, '-','color', [0.6 0.6, 0.6], LineWidth=4)
ylabel('Mass-Univariate Amplitude [µV]')
xlim([-1 0.2])
ylim([-2 3])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
subplot(7,1,4)
hold on
histogram(-sRTms(f5),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[181/255, 0, 0], 'LineWidth', 0.5)
histogram(-sRTms(f1),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[107/255, 0, 0], 'LineWidth',0.5)
xlim([-1000 200])
set(gca,'visible','off')
subplot(7,1,[5:7])
hold on
plot(tps,mean(All_y_Pisauro_deconv(:,:,1),1) , '-','color', [107/255 0, 0], LineWidth=2)
plot(tps, mean(All_y_Pisauro_deconv(:,:,2),1), '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, sig_Pisauro_deconv*-1.5, '-','color', [0.6 0.6, 0.6], LineWidth=4);
ylabel('Unfold-Deconvolution Amplitude [µV]')
xlim([-1 0.2])
ylim([-2 3])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
set(gcf,'position',figsize)
export_fig(sprintf('%sPisauro_Reanalysis_sANDr', figpath),'-pdf','-painters');

%% slocked
Pisauro_s_time = Pisauro.time(contains(Pisauro.rt_split, 'slower=') & contains(Pisauro.type,'deconv') & contains(Pisauro.basisname, 'stimulus') & Pisauro.time> -0.2 & Pisauro.time< 1);
All_y_Pisauro_MASS = nan(length(unique(Pisauro_r_subs)),length(unique(Pisauro_s_time)), 2);
All_y_Pisauro_deconv = nan(length(unique(Pisauro_r_subs)),length(unique(Pisauro_s_time)), 2);

tps = unique(Pisauro_s_time);
for tp = 1:length(tps)
    All_Betas_Pisauro_MASS(:,tp) = Pisauro_s_MASS_beta(Pisauro_s_time==tps(tp));
    All_Betas_Pisauro_deconv(:,tp) = Pisauro_s_deconv_beta(Pisauro_s_time==tps(tp));
    All_y_Pisauro_MASS(:,tp, 1) = Pisauro_s_MASS.yhat(contains(Pisauro_s_MASS.rt_split, 'slower=')& Pisauro_s_MASS.time==tps(tp));
    All_y_Pisauro_MASS(:,tp, 2) = Pisauro_s_MASS.yhat(contains(Pisauro_s_MASS.rt_split, 'faster')& Pisauro_s_MASS.time==tps(tp));
    All_y_Pisauro_deconv(:,tp, 1) = Pisauro_s_deconv.yhat(contains(Pisauro_s_deconv.rt_split, 'slower=')& Pisauro_s_deconv.time==tps(tp));
    All_y_Pisauro_deconv(:,tp, 2) = Pisauro_s_deconv.yhat(contains(Pisauro_s_deconv.rt_split, 'faster')& Pisauro_s_deconv.time==tps(tp));
end


%% now do t-test
thresh = 0.05;
[sig_Pisauro_MASS, pM, ~, dM]=ttest(All_Betas_Pisauro_MASS, [], thresh, 'both');
All_t_Pisauro_MASS = dM.tstat;
[sig_Pisauro_deconv, pD, ~, dD]=ttest(All_Betas_Pisauro_deconv, [], thresh, 'both');
All_t_Pisauro_deconv = dD.tstat;

sig_Pisauro_MASS(sig_Pisauro_MASS == 0) = nan;
sig_Pisauro_deconv(sig_Pisauro_deconv == 0) = nan;

sigClusts_Pisauro= [tps(sig_Pisauro_deconv==1), pD(sig_Pisauro_deconv==1)'];

%% permutation correction
timeLabels = tps;
PermStats_Pisauro_MASS = clustercorrection(All_Betas_Pisauro_MASS,[1], thresh, numPerm,timeLabels);
sig_Pisauro_MASS = sum(PermStats_Pisauro_MASS.maps,3);
sig_Pisauro_MASS(sig_Pisauro_MASS == 0) = nan;
PermStats_Pisauro_Deconv = clustercorrection(All_Betas_Pisauro_deconv,[1], thresh, numPerm,timeLabels);
try
sig_Pisauro_deconv = sum(PermStats_Pisauro_Deconv.maps,3);
sig_Pisauro_deconv(sig_Pisauro_deconv == 0) = nan;
catch
sig_Pisauro_deconv = nan(size(sig_Pisauro_MASS));
end

%%
figure()
set(gcf,'Color', [1 1 1])
subplot(7,1,[1:3])
hold on
plot(tps,mean(All_y_Pisauro_MASS(:,:,1),1) , '-','color', [107/255 0, 0], LineWidth=2)
plot(tps, mean(All_y_Pisauro_MASS(:,:,2),1), '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, sig_Pisauro_MASS*-1.5, '-','color', [0.6 0.6, 0.6], LineWidth=4)
ylabel('Mass-Univariate Amplitude [µV]')
xlim([-0.2 1])
ylim([-2 3])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
subplot(7,1,4)
hold on
histogram(sRTms(f5),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[181/255, 0, 0], 'LineWidth', 0.5)
histogram(sRTms(f1),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[107/255, 0, 0], 'LineWidth',0.5)
xlim([-200 1000])
set(gca,'visible','off')
subplot(7,1,[5:7])
hold on
plot(tps,mean(All_y_Pisauro_deconv(:,:,1),1) , '-','color', [107/255 0, 0], LineWidth=2)
plot(tps, mean(All_y_Pisauro_deconv(:,:,2),1), '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, sig_Pisauro_deconv*-1.5, '-','color', [0.6 0.6, 0.6], LineWidth=4);
ylabel('Unfold-Deconvolution Amplitude [µV]')
xlim([-0.2 1])
ylim([-2 3])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
set(gcf,'position',figsize)
export_fig(sprintf('%sPisauro_Reanalysis_slocked_sANDr', figpath),'-pdf','-painters');

%% UCAP

% what I need: 
% behavioral data for computing yhats
load('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Analyses/Matlab/UCAP/Second_level/Export/behavTable.mat');
[~,idx] = unique(behavTable.VPNummer,'first');
vps = behavTable.VPNummer(sort(idx));
behavTable.Quality = zeros(size(behavTable.VPNummer));
behavTable.Quality(strcmp(behavTable.n_b,'normal')==1)=1;
behavTable.DevP = zeros(size(behavTable.VPNummer));
behavTable.DevP(strcmp(behavTable.DeviantPosRL,'re')==1)=1;
behavTablecor = behavTable(behavTable.accuracy==1,:);

% filters for generating averages
intactleft = find(behavTablecor.Quality==1 & behavTablecor.DevP==0);
intactright = find(behavTablecor.Quality==1 & behavTablecor.DevP==1);

blurrleft = find(behavTablecor.Quality==0 & behavTablecor.DevP==0 );
blurrright = find(behavTablecor.Quality==0 & behavTablecor.DevP==1 );

% t-values - grab from 

% Unfold
load('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Analyses/Matlab/UCAP/Second_level/Export/regcoeffs.mat', 'sbasicCoeffs', 'rbasicCoeffs', 'chanlocs');
% y-hats for 
Pz = find(strcmp({chanlocs(:).labels},'Pz'));
tps = rbasicCoeffs.times;
% I'm grabbing the quality effect, because that's where the main difference
% is
All_Betas_UCAP_deconv = squeeze(rbasicCoeffs.data(Pz, :, 3, :))';
Deconv_r_betas = rbasicCoeffs;

clear yHat 
clear allYHats
allYHats = [];
lastn = 0;
curbetas = Deconv_r_betas.data; % let's do one first - first Valuation component 
for subj = 1:length(vps)
s_id = vps(subj);
clear yHat
        bSubData = behavTable(behavTable.VPNummer==s_id & behavTable.accuracy~= 0,:); % get a subjects data, I think I already removed non-responses, but whatever
        xMat=meanCentX([ones(size(bSubData,1),1),bSubData.DevP , bSubData.Quality, bSubData.DevP.*bSubData.Quality]); %zeros(size(bSubData,1),1) now I'm setting the intercept and the variable of interest to zero
        for t = 1: size(curbetas,2)    
        beta = squeeze(curbetas(Pz,t,:, subj)); 
        yHat(t, :)= [xMat]*beta;    
        end
        allYHats(:,lastn+1:lastn+size(yHat,2)) = yHat; 
        lastn = size(allYHats,2);
end

%% compute means
ERPintactleft = nanmean(allYHats(:,intactleft ),2);
ERPintactright = nanmean(allYHats(:,intactright ),2);
ERPblurrleft = nanmean(allYHats(:,blurrleft ),2);
ERPblurrright = nanmean(allYHats(:,blurrright ),2);%allYHats;

All_y_UCAP_deconv = [ERPintactleft'; ERPintactright'; ERPblurrleft'; ERPblurrright'];

%All_y_UCAP_deconv =
% MASS
load('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Analyses/Matlab/UCAP/Second_level/Export/regcoeffs_no_dc.mat')
All_Betas_UCAP_MASS = squeeze(rbasicCoeffs.data(Pz, :, 3, :))';
MASS_r_betas = rbasicCoeffs;


clear yHat 
clear allYHats
allYHats = [];
lastn = 0;
curbetas = MASS_r_betas.data; % let's do one first - first Valuation component 
for subj = 1:length(vps)
s_id = vps(subj);
clear yHat
        bSubData = behavTable(behavTable.VPNummer==s_id & behavTable.accuracy~= 0,:); % get a subjects data, I think I already removed non-responses, but whatever
        xMat=meanCentX([ones(size(bSubData,1),1),bSubData.DevP , bSubData.Quality, bSubData.DevP.*bSubData.Quality]); %zeros(size(bSubData,1),1) now I'm setting the intercept and the variable of interest to zero
        for t = 1: size(curbetas,2)    
        beta = squeeze(curbetas(Pz,t,:, subj)); 
        yHat(t, :)= [xMat]*beta;    
        end
        allYHats(:,lastn+1:lastn+size(yHat,2)) = yHat; 
        lastn = size(allYHats,2);
end

%% compute means
ERPintactleft = nanmean(allYHats(:,intactleft ),2);
ERPintactright = nanmean(allYHats(:,intactright ),2);
ERPblurrleft = nanmean(allYHats(:,blurrleft ),2);
ERPblurrright = nanmean(allYHats(:,blurrright ),2);%allYHats;

All_y_UCAP_MASS = [ERPintactleft'; ERPintactright'; ERPblurrleft'; ERPblurrright'];

%%
thresh = 0.05;
[sig_UCAP_MASS, pM, ~, dM]=ttest(All_Betas_UCAP_MASS, [], thresh, 'both');
All_t_UCAP_MASS = dM.tstat;
[sig_UCAP_deconv, pD, ~, dD]=ttest(All_Betas_UCAP_deconv, [], thresh, 'both');
All_t_UCAP_deconv = dD.tstat;

sig_UCAP_MASS(sig_UCAP_MASS == 0) = nan;
sig_UCAP_deconv(sig_UCAP_deconv == 0) = nan;

%% permutation correction
timeLabels = tps;
PermStats_UCAP_MASS = clustercorrection(All_Betas_UCAP_MASS,[1], thresh, numPerm,timeLabels);
sig_UCAP_MASS = sum(PermStats_UCAP_MASS.maps,3);
sig_UCAP_MASS(sig_UCAP_MASS == 0) = nan;
PermStats_UCAP_Deconv = clustercorrection(All_Betas_UCAP_deconv,[1], thresh, numPerm,timeLabels);
try
sig_UCAP_deconv = sum(PermStats_UCAP_Deconv.maps,3);
sig_UCAP_deconv(sig_UCAP_deconv == 0) = nan;
catch
sig_UCAP_deconv = nan(size(sig_UCAP_MASS));
end



%%
figure()
set(gcf,'Color', [1 1 1])
hold on
Leg1 = plot(tps, All_t_UCAP_MASS, 'k-', LineWidth=2);
plot(tps, sig_UCAP_MASS*-8, 'k-', LineWidth=4)
Leg2 = plot(tps, All_t_UCAP_deconv, '-','color', [0.6 0.6, 0.6], LineWidth=2);
plot(tps, sig_UCAP_deconv*-7.5, '-','color', [0.6 0.6, 0.6], LineWidth=4);
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'Mass-Univariate', 'Unfold-Deconvolution');% , Leg3, Leg4
set(gca,'TickDir','out', 'Box', 'off')
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14, 'position',[.25,.25,lp(3:4)])%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')

%% let's just plot the conditions:

figure()
set(gcf,'Color', [1 1 1])
hold on
subplot(7,1,[1:3])
hold on
Leg1=plot(tps,All_y_UCAP_MASS(1,:) , '-','color', [181/255, 0, 0], LineWidth=2);
Leg3=plot(tps, All_y_UCAP_MASS(3,:), '-','color', [107/255, 0, 0], LineWidth=2);
Leg2=plot(tps,All_y_UCAP_MASS(2,:) , '-.','color', [181/255, 0, 0], LineWidth=2);
Leg4=plot(tps, All_y_UCAP_MASS(4,:), '-.','color', [107/255, 0, 0], LineWidth=2);
plot(tps, sig_UCAP_MASS*-1.7, '-','color', [0.6 0.6, 0.6], LineWidth=4)
ylabel('Mass-Univariate Amplitude [µV]')
xlim([-1 0.2])
ylim([-2 5])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
h = legend([Leg1, Leg2, Leg3, Leg4],'intact left', 'intact right', 'blurr left', 'blurr right');% , Leg3, Leg4
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 12, 'position',[.25,.75,lp(3:4)])%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')
subplot(7,1,4)
hold on
histogram(-behavTablecor.RT(intactleft),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[181/255, 0, 0], 'LineWidth', 0.5)
histogram(-behavTablecor.RT(intactright),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[181/255, 0, 0], 'LineStyle', '-.', 'LineWidth', 0.5)
histogram(-behavTablecor.RT(blurrleft),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[107/255, 0, 0], 'LineWidth',0.5)
histogram(-behavTablecor.RT(blurrright),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[107/255, 0, 0],'LineStyle', '-.', 'LineWidth', 0.5)
xlim([-1000 200])
set(gca,'visible','off')
subplot(7,1,[5:7])
hold on
plot(tps,All_y_UCAP_deconv(1,:) , '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, All_y_UCAP_deconv(3,:), '-','color', [107/255 0, 0], LineWidth=2)
plot(tps,All_y_UCAP_deconv(2,:) , '-.','color', [181/255 0, 0], LineWidth=2)
plot(tps, All_y_UCAP_deconv(4,:), '-.','color', [107/255 0, 0], LineWidth=2)
plot(tps, sig_UCAP_deconv*-1.7, '-','color', [0.6 0.6, 0.6], LineWidth=4);
ylabel('Unfold-Deconvolution Amplitude [µV]')
xlim([-1 0.2])
ylim([-2 5])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
set(gcf,'position',figsize)
export_fig(sprintf('%sUCAP_Reanalysis_sANDr', figpath),'-pdf','-painters');

%% slocked
load('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Analyses/Matlab/UCAP/Second_level/Export/regcoeffs.mat', 'sbasicCoeffs', 'rbasicCoeffs', 'chanlocs');

tps = sbasicCoeffs.times;
% I'm grabbing the quality effect, because that's where the main difference
% is
All_Betas_UCAP_deconv = squeeze(sbasicCoeffs.data(Pz, :, 3, :))';
Deconv_s_betas = sbasicCoeffs;

clear yHat 
clear allYHats
allYHats = [];
lastn = 0;
curbetas = Deconv_s_betas.data; % let's do one first - first Valuation component 
for subj = 1:length(vps)
s_id = vps(subj);
clear yHat
        bSubData = behavTable(behavTable.VPNummer==s_id & behavTable.accuracy~= 0,:); % get a subjects data, I think I already removed non-responses, but whatever
        xMat=meanCentX([ones(size(bSubData,1),1),bSubData.DevP , bSubData.Quality, bSubData.DevP.*bSubData.Quality]); %zeros(size(bSubData,1),1) now I'm setting the intercept and the variable of interest to zero
        for t = 1: size(curbetas,2)    
        beta = squeeze(curbetas(Pz,t,:, subj)); 
        yHat(t, :)= [xMat]*beta;    
        end
        allYHats(:,lastn+1:lastn+size(yHat,2)) = yHat; 
        lastn = size(allYHats,2);
end

%% compute means
ERPintactleft = nanmean(allYHats(:,intactleft ),2);
ERPintactright = nanmean(allYHats(:,intactright ),2);
ERPblurrleft = nanmean(allYHats(:,blurrleft ),2);
ERPblurrright = nanmean(allYHats(:,blurrright ),2);%allYHats;

All_y_UCAP_deconv = [ERPintactleft'; ERPintactright'; ERPblurrleft'; ERPblurrright'];


%%
load('~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Analyses/Matlab/UCAP/Second_level/Export/regcoeffs_no_dc.mat')
All_Betas_UCAP_MASS = squeeze(sbasicCoeffs.data(Pz, :, 3, :))';
MASS_s_betas = sbasicCoeffs;


clear yHat 
clear allYHats
allYHats = [];
lastn = 0;
curbetas = MASS_s_betas.data; % let's do one first - first Valuation component 
for subj = 1:length(vps)
s_id = vps(subj);
clear yHat
        bSubData = behavTable(behavTable.VPNummer==s_id & behavTable.accuracy~= 0,:); % get a subjects data, I think I already removed non-responses, but whatever
        xMat=meanCentX([ones(size(bSubData,1),1),bSubData.DevP , bSubData.Quality, bSubData.DevP.*bSubData.Quality]); %zeros(size(bSubData,1),1) now I'm setting the intercept and the variable of interest to zero
        for t = 1: size(curbetas,2)    
        beta = squeeze(curbetas(Pz,t,:, subj)); 
        yHat(t, :)= [xMat]*beta;    
        end
        allYHats(:,lastn+1:lastn+size(yHat,2)) = yHat; 
        lastn = size(allYHats,2);
end

%% compute means
ERPintactleft = nanmean(allYHats(:,intactleft ),2);
ERPintactright = nanmean(allYHats(:,intactright ),2);
ERPblurrleft = nanmean(allYHats(:,blurrleft ),2);
ERPblurrright = nanmean(allYHats(:,blurrright ),2);%allYHats;

All_y_UCAP_MASS = [ERPintactleft'; ERPintactright'; ERPblurrleft'; ERPblurrright'];

%%
thresh = 0.05;
[sig_UCAP_MASS, pM, ~, dM]=ttest(All_Betas_UCAP_MASS, [], thresh, 'both');
All_t_UCAP_MASS = dM.tstat;
[sig_UCAP_deconv, pD, ~, dD]=ttest(All_Betas_UCAP_deconv, [], thresh, 'both');
All_t_UCAP_deconv = dD.tstat;

sig_UCAP_MASS(sig_UCAP_MASS == 0) = nan;
sig_UCAP_deconv(sig_UCAP_deconv == 0) = nan;

%% permutation correction
timeLabels = tps;
PermStats_UCAP_MASS = clustercorrection(All_Betas_UCAP_MASS,[1], thresh, numPerm,timeLabels);
sig_UCAP_MASS = sum(PermStats_UCAP_MASS.maps,3);
sig_UCAP_MASS(sig_UCAP_MASS == 0) = nan;
PermStats_UCAP_Deconv = clustercorrection(All_Betas_UCAP_deconv,[1], thresh, numPerm,timeLabels);
try
sig_UCAP_deconv = sum(PermStats_UCAP_Deconv.maps,3);
sig_UCAP_deconv(sig_UCAP_deconv == 0) = nan;
catch
sig_UCAP_deconv = nan(size(sig_UCAP_MASS));
end



%%
figure()
set(gcf,'Color', [1 1 1])
hold on
subplot(7,1,[1:3])
hold on
Leg1=plot(tps,All_y_UCAP_MASS(1,:) , '-','color', [181/255, 0, 0], LineWidth=2);
Leg3=plot(tps, All_y_UCAP_MASS(3,:), '-','color', [107/255, 0, 0], LineWidth=2);
Leg2=plot(tps,All_y_UCAP_MASS(2,:) , '-.','color', [181/255, 0, 0], LineWidth=2);
Leg4=plot(tps, All_y_UCAP_MASS(4,:), '-.','color', [107/255, 0, 0], LineWidth=2);
plot(tps, sig_UCAP_MASS*-1.7, '-','color', [0.6 0.6, 0.6], LineWidth=4)
ylabel('Mass-Univariate Amplitude [µV]')
xlim([-0.2 1])
ylim([-3 5])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
h = legend([Leg1, Leg2, Leg3, Leg4],'intact left', 'intact right', 'blurr left', 'blurr right');% , Leg3, Leg4
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 12, 'position',[.25,.75,lp(3:4)])%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')
subplot(7,1,4)
hold on
histogram(behavTablecor.RT(intactleft),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[181/255, 0, 0], 'LineWidth', 0.5)
histogram(behavTablecor.RT(intactright),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[181/255, 0, 0], 'LineStyle', '-.', 'LineWidth', 0.5)
histogram(behavTablecor.RT(blurrleft),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[107/255, 0, 0], 'LineWidth',0.5)
histogram(behavTablecor.RT(blurrright),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[107/255, 0, 0],'LineStyle', '-.', 'LineWidth', 0.5)
xlim([-200 1000])
set(gca,'visible','off')
subplot(7,1,[5:7])
hold on
plot(tps,All_y_UCAP_deconv(1,:) , '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, All_y_UCAP_deconv(3,:), '-','color', [107/255 0, 0], LineWidth=2)
plot(tps,All_y_UCAP_deconv(2,:) , '-.','color', [181/255 0, 0], LineWidth=2)
plot(tps, All_y_UCAP_deconv(4,:), '-.','color', [107/255 0, 0], LineWidth=2)
plot(tps, sig_UCAP_deconv*-1.7, '-','color', [0.6 0.6, 0.6], LineWidth=4);
ylabel('Unfold-Deconvolution Amplitude [µV]')
xlim([-0.2 1])
ylim([-3 5])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
set(gcf,'position',figsize)
export_fig(sprintf('%sUCAP_Reanalysis_slocked_sANDr', figpath),'-pdf','-painters');



%% Same for Steinemann
Steinemann = readtable(sprintf('%sSteinemann_2018/preprocessing_for_Unfold/post_unfold/steinemann/2023-11-7_2023_romy_steinemann_filteredrtsplit-both', basepath));
% temporarily deal with the saving issue:
Steinemann = Steinemann(contains(Steinemann.mean, 'high') & contains(Steinemann.var, 'high'), :);

load(sprintf('%sSteinemann_2018/preprocessing_for_Unfold/analyses/AllSubRTs.mat', basepath));
% separate out: stimulus vs response; deconv vs epoched

Steinemann_s_deconv = Steinemann(contains(Steinemann.type,'deconv') & contains(Steinemann.basisname, 'Stimulus') & Steinemann.time> -0.2 & Steinemann.time< 1,:);
Steinemann_r_deconv = Steinemann(contains(Steinemann.type,'deconv') & contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time< 0.2,:);

Steinemann_s_MASS = Steinemann(contains(Steinemann.type,'epoched') & contains(Steinemann.basisname, 'Stimulus') & Steinemann.time> -0.2 & Steinemann.time< 1,:);
Steinemann_r_MASS = Steinemann(contains(Steinemann.type,'epoched') & contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time< 0.2,:);

Steinemann_r_deconv_beta = Steinemann.yhat(contains(Steinemann.rt_split, 'slower=') & contains(Steinemann.type,'deconv') & contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time< 0.2)...
    - Steinemann.yhat(contains(Steinemann.rt_split, 'faster') & contains(Steinemann.type,'deconv') & contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time<0.2);

Steinemann_r_MASS_beta = Steinemann.yhat(contains(Steinemann.rt_split, 'slower=') & contains(Steinemann.type,'epoched') & contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time< 0.2)...
    - Steinemann.yhat(contains(Steinemann.rt_split, 'faster') & contains(Steinemann.type,'epoched') & contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time< 0.2);

Steinemann_s_deconv_beta = Steinemann.yhat(contains(Steinemann.rt_split, 'slower=') & contains(Steinemann.type,'deconv') & contains(Steinemann.basisname, 'Stimulus') & Steinemann.time> -0.2 & Steinemann.time< 1)...
    - Steinemann.yhat(contains(Steinemann.rt_split, 'faster') & contains(Steinemann.type,'deconv') & contains(Steinemann.basisname, 'Stimulus') & Steinemann.time> -0.2 & Steinemann.time< 1);

Steinemann_s_MASS_beta = Steinemann.yhat(contains(Steinemann.rt_split, 'slower=') & contains(Steinemann.type,'epoched') & contains(Steinemann.basisname, 'Stimulus') & Steinemann.time> -0.2 & Steinemann.time< 1)...
    - Steinemann.yhat(contains(Steinemann.rt_split, 'faster') & contains(Steinemann.type,'epoched') & contains(Steinemann.basisname, 'Stimulus') & Steinemann.time> -0.2 & Steinemann.time< 1);


Steinemann_r_time = Steinemann.time(contains(Steinemann.rt_split, 'slower=') & contains(Steinemann.type,'deconv') & contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time< 0.2);

Steinemann_r_subs = Steinemann.subject(contains(Steinemann.rt_split, 'slower=') & contains(Steinemann.type,'deconv') & contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time< 0.2);

%% goal: compute t-test at all time-points
% sub-goal: make matrix where columns are Subjects and timepoints are
% columns

All_Betas_Steinemann_MASS = nan(length(unique(Steinemann_r_subs)),length(unique(Steinemann_r_time)));
All_Betas_Steinemann_deconv = nan(length(unique(Steinemann_r_subs)),length(unique(Steinemann_r_time)));

All_y_Steinemann_MASS = nan(length(unique(Steinemann_r_subs)),length(unique(Steinemann_r_time)), 2);
All_y_Steinemann_deconv = nan(length(unique(Steinemann_r_subs)),length(unique(Steinemann_r_time)), 2);

tps = unique(Steinemann_r_time);
for tp = 1:length(tps)
    All_Betas_Steinemann_MASS(:,tp) = Steinemann_r_MASS_beta(Steinemann_r_time==tps(tp));
    All_Betas_Steinemann_deconv(:,tp) = Steinemann_r_deconv_beta(Steinemann_r_time==tps(tp));
    All_y_Steinemann_MASS(:,tp, 1) = Steinemann_r_MASS.yhat(contains(Steinemann_r_MASS.rt_split, 'slower=')& Steinemann_r_MASS.time==tps(tp));
    All_y_Steinemann_MASS(:,tp, 2) = Steinemann_r_MASS.yhat(contains(Steinemann_r_MASS.rt_split, 'faster')& Steinemann_r_MASS.time==tps(tp));
    All_y_Steinemann_deconv(:,tp, 1) = Steinemann_r_deconv.yhat(contains(Steinemann_r_deconv.rt_split, 'slower=')& Steinemann_r_deconv.time==tps(tp));
    All_y_Steinemann_deconv(:,tp, 2) = Steinemann_r_deconv.yhat(contains(Steinemann_r_deconv.rt_split, 'faster')& Steinemann_r_deconv.time==tps(tp));
end

%% now do t-test
thresh = 0.05;
[sig_Steinemann_MASS, pM, ~, dM]=ttest(All_Betas_Steinemann_MASS, [], thresh, 'both');
All_t_Steinemann_MASS = dM.tstat;
[sig_Steinemann_deconv, pD, ~, dD]=ttest(All_Betas_Steinemann_deconv, [], thresh, 'both');
All_t_Steinemann_deconv = dD.tstat;

sig_Steinemann_MASS(sig_Steinemann_MASS == 0) = nan;
sig_Steinemann_deconv(sig_Steinemann_deconv == 0) = nan;

sigClusts_Steinemann= [tps(sig_Steinemann_deconv==1), pD(sig_Steinemann_deconv==1)'];

%% permutation correction
timeLabels = tps;
PermStats_Steinemann_MASS = clustercorrection(All_Betas_Steinemann_MASS,[1], thresh, numPerm,timeLabels);
sig_Steinemann_MASS = sum(PermStats_Steinemann_MASS.maps,3);
sig_Steinemann_MASS(sig_Steinemann_MASS == 0) = nan;
PermStats_Steinemann_Deconv = clustercorrection(All_Betas_Steinemann_deconv,[1], thresh, numPerm,timeLabels);
try
sig_Steinemann_deconv = sum(PermStats_Steinemann_Deconv.maps,3);
sig_Steinemann_deconv(sig_Steinemann_deconv == 0) = nan;
catch
sig_Steinemann_deconv = nan(size(sig_Steinemann_MASS));
end

%%
figure()
set(gcf,'Color', [1 1 1])
hold on
Leg1 = plot(tps, All_t_Steinemann_MASS, 'k-', LineWidth=2);
plot(tps, sig_Steinemann_MASS*-8, 'k-', LineWidth=4)
Leg2 = plot(tps, All_t_Steinemann_deconv, '-','color', [0.6 0.6, 0.6], LineWidth=2);
plot(tps, sig_Steinemann_deconv*-7.5, '-','color', [0.6 0.6, 0.6], LineWidth=4);
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'Mass-Univariate', 'Unfold-Deconvolution');% , Leg3, Leg4
set(gca,'TickDir','out', 'Box', 'off')
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14, 'position',[.25,.25,lp(3:4)])%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')

%% let's just plot the conditions:

figure()
set(gcf,'Color', [1 1 1])
subplot(7,1,[1:3])
hold on
plot(tps,mean(All_y_Steinemann_MASS(:,:,1),1) , '-','color', [107/255 0, 0], LineWidth=2)
plot(tps, mean(All_y_Steinemann_MASS(:,:,2),1), '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, sig_Steinemann_MASS*-1.5, '-','color', [0.6 0.6, 0.6], LineWidth=4)
ylabel('Mass-Univariate Amplitude [µV]')
xlim([-1 0.2])
ylim([-2 5])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
subplot(7,1,4)
hold on
histogram(-allSubRTs(allSubRTs < median(allSubRTs))*1000,'BinWidth',50, 'FaceColor','none', 'EdgeColor',[181/255, 0, 0], 'LineWidth', 0.5)
histogram(-allSubRTs(allSubRTs > median(allSubRTs))*1000,'BinWidth',50, 'FaceColor','none', 'EdgeColor',[107/255, 0, 0], 'LineWidth',0.5)
xlim([-1000 200])
set(gca,'visible','off')
subplot(7,1,[5:7])
hold on
plot(tps,mean(All_y_Steinemann_deconv(:,:,1),1) , '-','color', [107/255 0, 0], LineWidth=2)
plot(tps, mean(All_y_Steinemann_deconv(:,:,2),1), '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, sig_Steinemann_deconv*-1.5, '-','color', [0.6 0.6, 0.6], LineWidth=4);
ylabel('Unfold-Deconvolution Amplitude [µV]')
xlim([-1 0.2])
ylim([-2 5])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
set(gcf,'position',figsize)
export_fig(sprintf('%sSteinemann_Reanalysis_sANDr', figpath),'-pdf','-painters');

%%

Steinemann_s_time = Steinemann.time(contains(Steinemann.rt_split, 'slower=') & contains(Steinemann.type,'deconv') & contains(Steinemann.basisname, 'Stimulus') & Steinemann.time> -0.2 & Steinemann.time< 1);
All_y_Steinemann_MASS = nan(length(unique(Steinemann_r_subs)),length(unique(Steinemann_s_time)), 2);
All_y_Steinemann_deconv = nan(length(unique(Steinemann_r_subs)),length(unique(Steinemann_s_time)), 2);

tps = unique(Steinemann_s_time);
for tp = 1:length(tps)
    All_Betas_Steinemann_MASS(:,tp) = Steinemann_s_MASS_beta(Steinemann_s_time==tps(tp));
    All_Betas_Steinemann_deconv(:,tp) = Steinemann_s_deconv_beta(Steinemann_s_time==tps(tp));
    All_y_Steinemann_MASS(:,tp, 1) = Steinemann_s_MASS.yhat(contains(Steinemann_s_MASS.rt_split, 'slower=')& Steinemann_s_MASS.time==tps(tp));
    All_y_Steinemann_MASS(:,tp, 2) = Steinemann_s_MASS.yhat(contains(Steinemann_s_MASS.rt_split, 'faster')& Steinemann_s_MASS.time==tps(tp));
    All_y_Steinemann_deconv(:,tp, 1) = Steinemann_s_deconv.yhat(contains(Steinemann_s_deconv.rt_split, 'slower=')& Steinemann_s_deconv.time==tps(tp));
    All_y_Steinemann_deconv(:,tp, 2) = Steinemann_s_deconv.yhat(contains(Steinemann_s_deconv.rt_split, 'faster')& Steinemann_s_deconv.time==tps(tp));
end

%% now do t-test
thresh = 0.05;
[sig_Steinemann_MASS, pM, ~, dM]=ttest(All_Betas_Steinemann_MASS, [], thresh, 'both');
All_t_Steinemann_MASS = dM.tstat;
[sig_Steinemann_deconv, pD, ~, dD]=ttest(All_Betas_Steinemann_deconv, [], thresh, 'both');
All_t_Steinemann_deconv = dD.tstat;

sig_Steinemann_MASS(sig_Steinemann_MASS == 0) = nan;
sig_Steinemann_deconv(sig_Steinemann_deconv == 0) = nan;

sigClusts_Steinemann= [tps(sig_Steinemann_deconv==1), pD(sig_Steinemann_deconv==1)'];

%% permutation correction
timeLabels = tps;
PermStats_Steinemann_MASS = clustercorrection(All_Betas_Steinemann_MASS,[1], thresh, numPerm,timeLabels);
sig_Steinemann_MASS = sum(PermStats_Steinemann_MASS.maps,3);
sig_Steinemann_MASS(sig_Steinemann_MASS == 0) = nan;
PermStats_Steinemann_Deconv = clustercorrection(All_Betas_Steinemann_deconv,[1], thresh, numPerm,timeLabels);
try
sig_Steinemann_deconv = sum(PermStats_Steinemann_Deconv.maps,3);
sig_Steinemann_deconv(sig_Steinemann_deconv == 0) = nan;
catch
sig_Steinemann_deconv = nan(size(sig_Steinemann_MASS));
end
%%
figure()
set(gcf,'Color', [1 1 1])
subplot(7,1,[1:3])
hold on
plot(tps,mean(All_y_Steinemann_MASS(:,:,1),1) , '-','color', [107/255 0, 0], LineWidth=2)
plot(tps, mean(All_y_Steinemann_MASS(:,:,2),1), '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, sig_Steinemann_MASS*-1.5, '-','color', [0.6 0.6, 0.6], LineWidth=4)
ylabel('Mass-Univariate Amplitude [µV]')
xlim([-0.2 1])
ylim([-2 5])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
subplot(7,1,4)
hold on
histogram(allSubRTs(allSubRTs < median(allSubRTs))*1000,'BinWidth',50, 'FaceColor','none', 'EdgeColor',[181/255, 0, 0], 'LineWidth', 0.5)
histogram(allSubRTs(allSubRTs > median(allSubRTs))*1000,'BinWidth',50, 'FaceColor','none', 'EdgeColor',[107/255, 0, 0], 'LineWidth',0.5)
xlim([-200 1000])
set(gca,'visible','off')
subplot(7,1,[5:7])
hold on
plot(tps,mean(All_y_Steinemann_deconv(:,:,1),1) , '-','color', [107/255 0, 0], LineWidth=2)
plot(tps, mean(All_y_Steinemann_deconv(:,:,2),1), '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, sig_Steinemann_deconv*-1.5, '-','color', [0.6 0.6, 0.6], LineWidth=4);
ylabel('Unfold-Deconvolution Amplitude [µV]')
xlim([-0.2 1])
ylim([-2 5])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
set(gcf,'position',figsize)
export_fig(sprintf('%sSteinemann_Reanalysis_slocked_sANDr', figpath),'-pdf','-painters');


%% Same for Boldt
Boldt = readtable(sprintf('%sBoldt_et_al_2019/prep_Boldt_unfold/post_unfold/boldt/2023-11-7_2023_romy_boldtrtsplit-both', basepath));

% temporarily deal with the saving issue:
Boldt = Boldt(contains(Boldt.mean, 'high') & contains(Boldt.var, 'high'), :);

load(sprintf('%sBoldt_et_al_2019/prep_Boldt_unfold/analyses/allSubRTs.mat', basepath));
allSubRTs(allSubRTs<100)= allSubRTs(allSubRTs<100)*1000;
allSubRTs(allSubRTs<0) = [];
% separate out: stimulus vs response; deconv vs epoched

Boldt_s_deconv = Boldt(contains(Boldt.type,'deconv') & contains(Boldt.basisname, 'Stimulus') & Boldt.time> -0.2 & Boldt.time< 1,:);
Boldt_r_deconv = Boldt(contains(Boldt.type,'deconv') & contains(Boldt.basisname, 'Response') & Boldt.time> -1 & Boldt.time< 0.2,:);

Boldt_s_MASS = Boldt(contains(Boldt.type,'epoched') & contains(Boldt.basisname, 'Stimulus') & Boldt.time> -0.2 & Boldt.time< 1,:);
Boldt_r_MASS = Boldt(contains(Boldt.type,'epoched') & contains(Boldt.basisname, 'Response') & Boldt.time> -1 & Boldt.time< 0.2,:);

Boldt_r_deconv_beta = Boldt.yhat(contains(Boldt.rt_split, 'slower=') & contains(Boldt.type,'deconv') & contains(Boldt.basisname, 'Response') & Boldt.time> -1 & Boldt.time< 0.2)...
    - Boldt.yhat(contains(Boldt.rt_split, 'faster') & contains(Boldt.type,'deconv') & contains(Boldt.basisname, 'Response') & Boldt.time> -1 & Boldt.time<0.2);

Boldt_r_MASS_beta = Boldt.yhat(contains(Boldt.rt_split, 'slower=') & contains(Boldt.type,'epoched') & contains(Boldt.basisname, 'Response') & Boldt.time> -1 & Boldt.time< 0.2)...
    - Boldt.yhat(contains(Boldt.rt_split, 'faster') & contains(Boldt.type,'epoched') & contains(Boldt.basisname, 'Response') & Boldt.time> -1 & Boldt.time< 0.2);

Boldt_s_deconv_beta = Boldt.yhat(contains(Boldt.rt_split, 'slower=') & contains(Boldt.type,'deconv') & contains(Boldt.basisname, 'Stimulus') & Boldt.time> -0.2 & Boldt.time< 1)...
    - Boldt.yhat(contains(Boldt.rt_split, 'faster') & contains(Boldt.type,'deconv') & contains(Boldt.basisname, 'Stimulus') & Boldt.time> -0.2 & Boldt.time< 1);

Boldt_s_MASS_beta = Boldt.yhat(contains(Boldt.rt_split, 'slower=') & contains(Boldt.type,'epoched') & contains(Boldt.basisname, 'Stimulus') & Boldt.time> -0.2 & Boldt.time< 1)...
    - Boldt.yhat(contains(Boldt.rt_split, 'faster') & contains(Boldt.type,'epoched') & contains(Boldt.basisname, 'Stimulus') & Boldt.time> -0.2 & Boldt.time< 1);


Boldt_r_time = Boldt.time(contains(Boldt.rt_split, 'slower=') & contains(Boldt.type,'deconv') & contains(Boldt.basisname, 'Response') & Boldt.time> -1 & Boldt.time< 0.2);

Boldt_r_subs = Boldt.subject(contains(Boldt.rt_split, 'slower=') & contains(Boldt.type,'deconv') & contains(Boldt.basisname, 'Response') & Boldt.time> -1 & Boldt.time< 0.2);

%% goal: compute t-test at all time-points
% sub-goal: make matrix where columns are Subjects and timepoints are
% columns

All_Betas_Boldt_MASS = nan(length(unique(Boldt_r_subs)),length(unique(Boldt_r_time)));
All_Betas_Boldt_deconv = nan(length(unique(Boldt_r_subs)),length(unique(Boldt_r_time)));

All_y_Boldt_MASS = nan(length(unique(Boldt_r_subs)),length(unique(Boldt_r_time)), 2);
All_y_Boldt_deconv = nan(length(unique(Boldt_r_subs)),length(unique(Boldt_r_time)), 2);

tps = unique(Boldt_r_time);
for tp = 1:length(tps)
    All_Betas_Boldt_MASS(:,tp) = Boldt_r_MASS_beta(Boldt_r_time==tps(tp));
    All_Betas_Boldt_deconv(:,tp) = Boldt_r_deconv_beta(Boldt_r_time==tps(tp));
    All_y_Boldt_MASS(:,tp, 1) = Boldt_r_MASS.yhat(contains(Boldt_r_MASS.rt_split, 'slower=')& Boldt_r_MASS.time==tps(tp));
    All_y_Boldt_MASS(:,tp, 2) = Boldt_r_MASS.yhat(contains(Boldt_r_MASS.rt_split, 'faster')& Boldt_r_MASS.time==tps(tp));
    All_y_Boldt_deconv(:,tp, 1) = Boldt_r_deconv.yhat(contains(Boldt_r_deconv.rt_split, 'slower=')& Boldt_r_deconv.time==tps(tp));
    All_y_Boldt_deconv(:,tp, 2) = Boldt_r_deconv.yhat(contains(Boldt_r_deconv.rt_split, 'faster')& Boldt_r_deconv.time==tps(tp));
end

%% now do t-test
thresh = 0.05;
[sig_Boldt_MASS, pM, ~, dM]=ttest(All_Betas_Boldt_MASS, [], thresh, 'both');
All_t_Boldt_MASS = dM.tstat;
[sig_Boldt_deconv, pD, ~, dD]=ttest(All_Betas_Boldt_deconv, [], thresh, 'both');
All_t_Boldt_deconv = dD.tstat;

sig_Boldt_MASS(sig_Boldt_MASS == 0) = nan;
sig_Boldt_deconv(sig_Boldt_deconv == 0) = nan;

sigClusts_Boldt= [tps(sig_Boldt_deconv==1), pD(sig_Boldt_deconv==1)'];

%% permutation correction
timeLabels = tps;
PermStats_Boldt_MASS = clustercorrection(All_Betas_Boldt_MASS,[1], thresh, numPerm,timeLabels);
sig_Boldt_MASS = sum(PermStats_Boldt_MASS.maps,3);
sig_Boldt_MASS(sig_Boldt_MASS == 0) = nan;
PermStats_Boldt_Deconv = clustercorrection(All_Betas_Boldt_deconv,[1], thresh, numPerm,timeLabels);
try
sig_Boldt_deconv = sum(PermStats_Boldt_Deconv.maps,3);
sig_Boldt_deconv(sig_Boldt_deconv == 0) = nan;
catch
sig_Boldt_deconv = nan(size(sig_Boldt_MASS));
end

%%
figure()
set(gcf,'Color', [1 1 1])
hold on
Leg1 = plot(tps, All_t_Boldt_MASS, 'k-', LineWidth=2);
plot(tps, sig_Boldt_MASS*-8, 'k-', LineWidth=4)
Leg2 = plot(tps, All_t_Boldt_deconv, '-','color', [0.6 0.6, 0.6], LineWidth=2);
plot(tps, sig_Boldt_deconv*-7.5, '-','color', [0.6 0.6, 0.6], LineWidth=4);
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'Mass-Univariate', 'Unfold-Deconvolution');% , Leg3, Leg4
set(gca,'TickDir','out', 'Box', 'off')
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14, 'position',[.25,.25,lp(3:4)])%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')

%% let's just plot the conditions:

figure()
set(gcf,'Color', [1 1 1])
subplot(7,1,[1:3])
hold on
plot(tps,mean(All_y_Boldt_MASS(:,:,1),1) , '-','color', [107/255 0, 0], LineWidth=2)
plot(tps, mean(All_y_Boldt_MASS(:,:,2),1), '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, sig_Boldt_MASS*-2, '-','color', [0.6 0.6, 0.6], LineWidth=4)
ylabel('Mass-Univariate Amplitude [µV]')
xlim([-1 0.2])
ylim([-3 8])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
subplot(7,1,4)
hold on
histogram(-allSubRTs(allSubRTs < median(allSubRTs)),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[181/255, 0, 0], 'LineWidth', 0.5)
histogram(-allSubRTs(allSubRTs > median(allSubRTs)),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[107/255, 0, 0], 'LineWidth',0.5)
xlim([-1000 200])
set(gca,'visible','off')
subplot(7,1,[5:7])
hold on
plot(tps,mean(All_y_Boldt_deconv(:,:,1),1) , '-','color', [107/255 0, 0], LineWidth=2)
plot(tps, mean(All_y_Boldt_deconv(:,:,2),1), '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, sig_Boldt_deconv*-2, '-','color', [0.6 0.6, 0.6], LineWidth=4);
ylabel('Unfold-Deconvolution Amplitude [µV]')
xlim([-1 0.2])
ylim([-3 8])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
set(gcf,'position',figsize)
export_fig(sprintf('%sBoldt_Reanalysis_sANDr', figpath),'-pdf','-painters');

%%

Boldt_s_time = Boldt.time(contains(Boldt.rt_split, 'slower=') & contains(Boldt.type,'deconv') & contains(Boldt.basisname, 'Stimulus') & Boldt.time> -0.2 & Boldt.time< 1);
All_y_Boldt_MASS = nan(length(unique(Boldt_r_subs)),length(unique(Boldt_s_time)), 2);
All_y_Boldt_deconv = nan(length(unique(Boldt_r_subs)),length(unique(Boldt_s_time)), 2);

tps = unique(Boldt_s_time);
for tp = 1:length(tps)
    All_Betas_Boldt_MASS(:,tp) = Boldt_s_MASS_beta(Boldt_s_time==tps(tp));
    All_Betas_Boldt_deconv(:,tp) = Boldt_s_deconv_beta(Boldt_s_time==tps(tp));
    All_y_Boldt_MASS(:,tp, 1) = Boldt_s_MASS.yhat(contains(Boldt_s_MASS.rt_split, 'slower=')& Boldt_s_MASS.time==tps(tp));
    All_y_Boldt_MASS(:,tp, 2) = Boldt_s_MASS.yhat(contains(Boldt_s_MASS.rt_split, 'faster')& Boldt_s_MASS.time==tps(tp));
    All_y_Boldt_deconv(:,tp, 1) = Boldt_s_deconv.yhat(contains(Boldt_s_deconv.rt_split, 'slower=')& Boldt_s_deconv.time==tps(tp));
    All_y_Boldt_deconv(:,tp, 2) = Boldt_s_deconv.yhat(contains(Boldt_s_deconv.rt_split, 'faster')& Boldt_s_deconv.time==tps(tp));
end

%% now do t-test
thresh = 0.05;
[sig_Boldt_MASS, pM, ~, dM]=ttest(All_Betas_Boldt_MASS, [], thresh, 'both');
All_t_Boldt_MASS = dM.tstat;
[sig_Boldt_deconv, pD, ~, dD]=ttest(All_Betas_Boldt_deconv, [], thresh, 'both');
All_t_Boldt_deconv = dD.tstat;

sig_Boldt_MASS(sig_Boldt_MASS == 0) = nan;
sig_Boldt_deconv(sig_Boldt_deconv == 0) = nan;

sigClusts_Boldt= [tps(sig_Boldt_deconv==1), pD(sig_Boldt_deconv==1)'];

%% permutation correction
timeLabels = tps;
PermStats_Boldt_MASS = clustercorrection(All_Betas_Boldt_MASS,[1], thresh, numPerm,timeLabels);
sig_Boldt_MASS = sum(PermStats_Boldt_MASS.maps,3);
sig_Boldt_MASS(sig_Boldt_MASS == 0) = nan;
PermStats_Boldt_Deconv = clustercorrection(All_Betas_Boldt_deconv,[1], thresh, numPerm,timeLabels);
try
sig_Boldt_deconv = sum(PermStats_Boldt_Deconv.maps,3);
sig_Boldt_deconv(sig_Boldt_deconv == 0) = nan;
catch
sig_Boldt_deconv = nan(size(sig_Boldt_MASS));
end
%%
figure()
set(gcf,'Color', [1 1 1])
subplot(7,1,[1:3])
hold on
plot(tps,mean(All_y_Boldt_MASS(:,:,1),1) , '-','color', [107/255 0, 0], LineWidth=2)
plot(tps, mean(All_y_Boldt_MASS(:,:,2),1), '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, sig_Boldt_MASS*-2, '-','color', [0.6 0.6, 0.6], LineWidth=4)
ylabel('Mass-Univariate Amplitude [µV]')
xlim([-0.2 1])
ylim([-2 9])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
subplot(7,1,4)
hold on
histogram(allSubRTs(allSubRTs < median(allSubRTs)),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[181/255, 0, 0], 'LineWidth', 0.5)
histogram(allSubRTs(allSubRTs > median(allSubRTs)),'BinWidth',50, 'FaceColor','none', 'EdgeColor',[107/255, 0, 0], 'LineWidth',0.5)
xlim([-200 1000])
set(gca,'visible','off')
subplot(7,1,[5:7])
hold on
plot(tps,mean(All_y_Boldt_deconv(:,:,1),1) , '-','color', [107/255 0, 0], LineWidth=2)
plot(tps, mean(All_y_Boldt_deconv(:,:,2),1), '-','color', [181/255 0, 0], LineWidth=2)
plot(tps, sig_Boldt_deconv*-2, '-','color', [0.6 0.6, 0.6], LineWidth=4);
ylabel('Unfold-Deconvolution Amplitude [µV]')
xlim([-0.2 1])
ylim([-2 9])
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gca,'TickDir','out', 'Box', 'off')
set(gcf,'position',figsize)
export_fig(sprintf('%sBoldt_Reanalysis_slocked_sANDr', figpath),'-pdf','-painters');


