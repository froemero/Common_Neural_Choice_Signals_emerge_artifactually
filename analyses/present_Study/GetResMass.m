
% Intent: get information saved out for each component to easily write up
% results

%%
PATH = '/Users/Romy/Dropbox (Brown)/ShenhavLab/EEG_ressources/Experiments/BASB_EEG/';
FIGPATH = 'Figures/';
addpath(genpath('~/Dropbox (Brown)/ShenhavLab/EEG_ressources/EEGfunctions/'));
load(strcat(PATH,'/Data/Export/TFchanlocs.mat'))

clc;
close all;
isROI=0;

delete(strcat(PATH, 'Writing/', 'MASSRESULTS'))
diary(strcat(PATH, 'Writing/', 'MASSRESULTS'))
nameFolds= { 'WeightedT_Liking_Anxiety_Confidence','WeightedT_Chosen_Unchosen','WeightedT_VD_OV','WeightedT_PCValuation_PCChoiceDifficulty'};    
    
figsize= [0 0 25 5];

threshVals= {'0.001', '0.005'};
MAPLIM= [-3 3];
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
warning('off');
ColorVals = {[88/255, 116/255, 152/255] ,[255/255, 216/255, 0/255]; [101/255, 153/255, 255/255],[ 255/255, 156/255, 0/255]; [0/255, 114/255, 178/255],[213/255, 94/255, 0/255]; [0/255, 153/255, 204/255],[204/255, 121/255, 167/255]; [102/255, 250/255, 255/255],[204/255, 121/255, 167/255]; [102/255, 100/255, 255/255],[255/255, 100/255, 0/255]  };
for nmodels = 1: length(nameFolds)
    fprintf(strcat('Results for:',nameFolds{nmodels},'\n\n') )
    for currthresh = 1: length(threshVals)
        %%
        if ~isROI
        FIGNAME = strcat(nameFolds{nmodels},'_', threshVals{currthresh}, '_effect_summary');
        
        figure();
        set(gcf,'Color', [1 1 1]);
        
        hold on;
        end
        colcount=0;
        load(strcat(PATH, 'Data/Export/',nameFolds{nmodels}, '_',threshVals{currthresh},'/ROI_data.mat'))
        % get fieldnames
        fns = fieldnames(allROIs);
        fprintf('Results at cluster forming threshold %.3f\n\n',allROIs.(fns{1}).clustThresh)
        
        sconds={};
        rconds={};
        for ii = 1:length(fns)
            
            
            if (~ isempty(allROIs.(fns{ii}).sign) & isempty(strfind(fns{ii}, 'intercept')))
                colcount = colcount+1;
                fprintf('%s\n',fns{ii}(13:end))
                fprintf('%s - locked\n',fns{ii}(1))
                if ~isROI
                if strcmp(fns{ii}(1), 's')
                    subplot(1,8,[1 2 ])
                    set(gca,'TickDir','out');
                    hold on
                    xlim([0 1000])
                    xlabel('time [ms]')
                    ylim([-1 60])
                    yticks([10 30 50])
                    yticklabels({'posterior','central','frontal'})
                    ylabel('electrode position')
                elseif strcmp(fns{ii}(1), 'r')
                    subplot(1,8,[5 6 ])
                    set(gca,'TickDir','out');
                    hold on
                    xlim([-1000 0])
                    xlabel('time [ms]')
                    ylim([-1 60])
                    yticks([10 30 50])
                    yticklabels({'posterior','central','frontal'})
                    ylabel('electrode position')
                end
                end
                % for a given predictor go through all significant effects
                for ncomps = 1: size(allROIs.(fns{ii}).sign,1)
                    % for positive components
                    if (allROIs.(fns{ii}).sign(ncomps)==1)
                        fprintf('Positive cluster at p = %.3f\n', allROIs.(fns{ii}).p(ncomps))
                        if strcmp(fns{ii}(1), 's')
                            if isempty(sconds)
                                sconds= cellstr(sprintf( '%s\n positive',fns{ii}(13:end)));
                            else
                                sconds(end+1)= cellstr(sprintf( '%s\n positive',fns{ii}(13:end)));
                            end
                            
                        elseif strcmp(fns{ii}(1), 'r') 
                            if isempty(rconds)
                                rconds= cellstr(sprintf( '%s\n positive',fns{ii}(13:end)));
                                if (~isempty(sconds) & strcmp(rconds{1}, sconds{1}))
                                 colcount = colcount-1;   
                                end
                            else
                                rconds(end+1)= cellstr(sprintf( '%s\n positive',fns{ii}(13:end)));
                                
                            end
                            
                            
                        end
                        allROIelecs=[];
                        y1 = 61 - allROIs.(fns{ii}).peakChannel(ncomps);
                        if allROIs.(fns{ii}).peakChannel(ncomps)==61
                            y1 =  61 - 30;
                        end
                        y2 = y1-1;%
                        for nROIelecs= 1: length(allROIs.(fns{ii}).roielecs{ncomps})
                            allROIelecs= strcat(allROIelecs, TFchanlocs(allROIs.(fns{ii}).roielecs{ncomps}(nROIelecs)).labels, ', ');
                        end
                        fprintf('Electrodes in cluster: %s\n',allROIelecs(1:end-1))
                        fprintf('Cluster time: %d ms : %d ms \n',allROIs.(fns{ii}).roitimes{ncomps}(1), allROIs.(fns{ii}).roitimes{ncomps}(end))
                        x1 = allROIs.(fns{ii}).roitimes{ncomps}(1);
                        x2=  allROIs.(fns{ii}).roitimes{ncomps}(end);
                        fprintf('Peak: %s at %d ms \n\n',TFchanlocs(allROIs.(fns{ii}). peakChannel(ncomps)).labels,allROIs.(fns{ii}).peakTime(ncomps))
                        if (~isROI & ~strcmp(fns{ii}(1), 'p'))
                        x = [x1, x2, x2, x1, x1];
                        y = [y1, y1, y2, y2, y1];
                        fill(x, y, [ColorVals{colcount, 2}], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                        % else for negative components
                        end
                    else
                        fprintf('Negative cluster at p = %.3f\n', allROIs.(fns{ii}).p(ncomps))
                        if strcmp(fns{ii}(1), 's')
                            if isempty(sconds)
                                sconds= cellstr(sprintf( '%s \n negative',fns{ii}(13:end)));
                            else
                                sconds(end+1)= cellstr(sprintf( '%s \n negative',fns{ii}(13:end)));
                            end
                            
                        elseif strcmp(fns{ii}(1), 'r') 
                            if isempty(rconds)
                                rconds= cellstr(sprintf( '%s \n negative',fns{ii}(13:end)));
                            else
                                rconds(end+1)= cellstr(sprintf( '%s \n negative',fns{ii}(13:end)));
                            end
                            
                        end
                        allROIelecs=[];
                        y1 = 61 - allROIs.(fns{ii}).peakChannel(ncomps);
                        if allROIs.(fns{ii}).peakChannel(ncomps)==61
                            y1 =  61 - 30;
                        end
                        y2 = y1-1;
                        for nROIelecs= 1: length(allROIs.(fns{ii}).roielecs{ncomps})
                            allROIelecs= strcat(allROIelecs, TFchanlocs(allROIs.(fns{ii}).roielecs{ncomps}(nROIelecs)).labels, ', ');
                        end
                        fprintf('Electrodes in cluster: %s\n',allROIelecs(1:end-1))
                        fprintf('Cluster time: %d ms : %d ms \n',allROIs.(fns{ii}).roitimes{ncomps}(1), allROIs.(fns{ii}).roitimes{ncomps}(end))
                        x1 = allROIs.(fns{ii}).roitimes{ncomps}(1);
                        x2=  allROIs.(fns{ii}).roitimes{ncomps}(end);
                        fprintf('Peak: %s at %d ms \n\n',TFchanlocs(allROIs.(fns{ii}).peakChannel(ncomps)).labels,allROIs.(fns{ii}).peakTime(ncomps))
                        if(~isROI & ~strcmp(fns{ii}(1), 'p'))
                        x = [x1, x2, x2, x1, x1];
                        y = [y1, y1, y2, y2, y1];
                        fill(x, y, [ColorVals{colcount, 1}], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                        end
                    end
                    %% save figure of current component
                    if (~isROI & ~strcmp(fns{ii}(1), 'p'))
                    figure();
                    topoplot(mean(allROIs.(fns{ii}).Tmaps(:,:, ncomps),2)*allROIs.(fns{ii}).sign(ncomps),TFchanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'colormap', redblue, 'conv', 'on');
                    xlim([-0.7 0.7])
                    ylim([-0.7 0.7])
                    export_fig(fullfile(strcat(PATH,FIGPATH,nameFolds{nmodels}, '_',fns{ii},'_',threshVals{currthresh},'_',num2str(ncomps))),'-pdf','-painters');
                    close gcf
                    end
                end
                
            end
            
        end
        if ~isROI
        try
            subplot(1,8,[1 2 ])
            ax = gca;
            ax.YTickLabelRotation = 90;
            %sconds
            lh=legend(sconds, 'Location', 'eastoutside','Box','off');
            axP = get(gca,'Position');
            lp=get(lh,'position');
            set(lh, 'position',[axP(1)+axP(3)+0.01 , lp(2),0.1, 0.1])
            subplot(1,8,[5 6 ])
            ax = gca;
            ax.YTickLabelRotation = 90;
            %rconds
            lh=legend(rconds, 'Location', 'eastoutside','Box','off');
            axP = get(gca,'Position');
            lp=get(lh,'position');
            set(lh,  'position',[axP(1)+axP(3)+0.01,lp(2),0.1, 0.1]);
        catch
        end
        end
        fprintf('\n\n')
        if~isROI 
        subtitle(sprintf( '%s %s',strrep(nameFolds{nmodels}, '_', ' '),threshVals{currthresh}));
        set(gcf,'units','centimeters','position',figsize)
        export_fig(sprintf('%s%s%s', PATH,FIGPATH,FIGNAME),'-pdf','-painters');
        end
    end
end

diary off





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
