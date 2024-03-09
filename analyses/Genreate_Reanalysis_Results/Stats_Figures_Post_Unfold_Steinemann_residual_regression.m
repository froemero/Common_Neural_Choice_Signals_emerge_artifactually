
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


%% Same for Steinemann
Steinemann = readtable(sprintf('%sSteinemann_2018/preprocessing_for_Unfold/post_unfold/steinemann/2024-03-08_massunivariateresiduals.csv', basepath));

% table has data per subject, so we will need to average first

% ogSteinemann = Steinemann;
% 
% subjs = unique(Steinemann.subject);
% 
% SteinemannY = nan(length(ogSteinemann.subject(ogSteinemann.subject==subjs(1))), 1, length(subjs));
% for ii= 1:length(subjs)
%     
%     SteinemannY(:,:,ii)= ogSteinemann.yhat(ogSteinemann.subject==subjs(ii), :);
% end
% 
% Steinemann = ogSteinemann(ogSteinemann.subject==subjs(1),:);
% Steinemann.yhat = mean(SteinemannY,3);

% separate out: stimulus vs response; deconv vs epoched

Steinemann_s_MASS = Steinemann(contains(Steinemann.basisname, 'Stimulus') & Steinemann.time> -0.2 & Steinemann.time< 1,:);
Steinemann_r_MASS = Steinemann(contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time< 0.2,:);

Steinemann_r_MASS_beta = Steinemann.yhat(contains(Steinemann.rt_split, 'slower=') &  contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time< 0.2)...
    - Steinemann.yhat(contains(Steinemann.rt_split, 'faster') & contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time< 0.2);

Steinemann_s_MASS_beta = Steinemann.yhat(contains(Steinemann.rt_split, 'slower=') &  contains(Steinemann.basisname, 'Stimulus') & Steinemann.time> -0.2 & Steinemann.time< 1)...
    - Steinemann.yhat(contains(Steinemann.rt_split, 'faster') & contains(Steinemann.basisname, 'Stimulus') & Steinemann.time> -0.2 & Steinemann.time< 1);


Steinemann_r_time = Steinemann.time(contains(Steinemann.rt_split, 'slower=') & contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time< 0.2);

Steinemann_r_subs = Steinemann.subject(contains(Steinemann.rt_split, 'slower=') & contains(Steinemann.basisname, 'Response') & Steinemann.time> -1 & Steinemann.time< 0.2);

%% goal: compute t-test at all time-points
% sub-goal: make matrix where columns are Subjects and timepoints are
% columns

All_Betas_Steinemann_MASS = nan(length(unique(Steinemann_r_subs)),length(unique(Steinemann_r_time)));
All_y_Steinemann_MASS = nan(length(unique(Steinemann_r_subs)),length(unique(Steinemann_r_time)), 2);

tps = unique(Steinemann_r_time);
for tp = 1:length(tps)
    All_Betas_Steinemann_MASS(:,tp) = Steinemann_r_MASS_beta(Steinemann_r_time==tps(tp));
    All_y_Steinemann_MASS(:,tp, 1) = Steinemann_r_MASS.yhat(contains(Steinemann_r_MASS.rt_split, 'slower=')& Steinemann_r_MASS.time==tps(tp));
    All_y_Steinemann_MASS(:,tp, 2) = Steinemann_r_MASS.yhat(contains(Steinemann_r_MASS.rt_split, 'faster')& Steinemann_r_MASS.time==tps(tp));
end

%% now do t-test
thresh = 0.05;
[sig_Steinemann_MASS, pM, ~, dM]=ttest(All_Betas_Steinemann_MASS, [], thresh, 'both');
All_t_Steinemann_MASS = dM.tstat;

sig_Steinemann_MASS(sig_Steinemann_MASS == 0) = nan;

%% permutation correction
timeLabels = tps;
PermStats_Steinemann_MASS = clustercorrection(All_Betas_Steinemann_MASS,[1], thresh, numPerm,timeLabels);
try
sig_Steinemann_MASS = sum(PermStats_Steinemann_MASS.maps,3);
sig_Steinemann_MASS(sig_Steinemann_MASS == 0) = nan;
catch
sig_Steinemann_MASS = nan(size(PermStats_Steinemann_MASS, 1), size(PermStats_Steinemann_MASS, 2));
end


%%
figure()
set(gcf,'Color', [1 1 1])
hold on
Leg1 = plot(tps, All_t_Steinemann_MASS, 'k-', LineWidth=2);
plot(tps, sig_Steinemann_MASS*-8, 'k-', LineWidth=4)
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1],'Mass-Univariate');% , Leg3, Leg4
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
set(gcf,'position',figsize)
export_fig(sprintf('%sSteinemann_ResidualRegressionr', figpath),'-pdf','-painters');

%%

Steinemann_s_time = Steinemann.time(contains(Steinemann.rt_split, 'slower=') & contains(Steinemann.basisname, 'Stimulus') & Steinemann.time> -0.2 & Steinemann.time< 1);
All_y_Steinemann_MASS = nan(length(unique(Steinemann_r_subs)),length(unique(Steinemann_s_time)), 2);

tps = unique(Steinemann_s_time);
for tp = 1:length(tps)
    All_Betas_Steinemann_MASS(:,tp) = Steinemann_s_MASS_beta(Steinemann_s_time==tps(tp));
    All_y_Steinemann_MASS(:,tp, 1) = Steinemann_s_MASS.yhat(contains(Steinemann_s_MASS.rt_split, 'slower=')& Steinemann_s_MASS.time==tps(tp));
    All_y_Steinemann_MASS(:,tp, 2) = Steinemann_s_MASS.yhat(contains(Steinemann_s_MASS.rt_split, 'faster')& Steinemann_s_MASS.time==tps(tp));
end

%% now do t-test
thresh = 0.05;
[sig_Steinemann_MASS, pM, ~, dM]=ttest(All_Betas_Steinemann_MASS, [], thresh, 'both');
All_t_Steinemann_MASS = dM.tstat;

sig_Steinemann_MASS(sig_Steinemann_MASS == 0) = nan;

%% permutation correction
timeLabels = tps;
PermStats_Steinemann_MASS = clustercorrection(All_Betas_Steinemann_MASS,[1], thresh, numPerm,timeLabels);
try
sig_Steinemann_MASS = sum(PermStats_Steinemann_MASS.maps,3);
sig_Steinemann_MASS(sig_Steinemann_MASS == 0) = nan;
catch
sig_Steinemann_MASS = nan(size(PermStats_Steinemann_MASS,1),size(PermStats_Steinemann_MASS,2));
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
set(gcf,'position',figsize)
export_fig(sprintf('%sSteinemann_Reanalysis_slocked_sANDr', figpath),'-pdf','-painters');

