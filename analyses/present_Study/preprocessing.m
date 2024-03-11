%% ATTENTION, stuff has been modified for no hiogh pass filter test

%% preparation
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

PATH = '/Users/Romy/Dropbox (Brown)/ShenhavLab/EEG_ressources/Experiments/BASB_EEG/'; %edit to your needs 

SOURCEFILES = dir(strcat(PATH, 'Data/raw/*.vhdr')); %all brain vision header files in the folder
SUBJECTS = 1:numel(SOURCEFILES);
%%
for s =1:numel(SOURCEFILES)%SUBJECTS
%% load Data
    
    fn   = SOURCEFILES(s).name(1:4);
    s_id = str2double(SOURCEFILES(s).name(3:4)) % reads in subject ID and converts it to double (does not work for '01')
    if isnan(s_id) % if operation returns NaN (because of zero)
    s_id = str2double(SOURCEFILES(s).name(4)) % then convert only second part of the string (2-9 in our Data)  
    end
    
    
    %%
    EEG = pop_loadbv(sprintf('%sData/raw/',PATH),sprintf('%s',SOURCEFILES(s).name));


%% do OC 

%add empty channel Cz
    EEG.data(end+1,:) = 0;
    EEG.nbchan = size(EEG.data,1);
    EEG.chanlocs(end+1).labels = 'Cz';
%%
  
    %load chanloc information
    %EEG = pop_chanedit(EEG, 'lookup',sprintf('%sOccular_Correction/werfen_BESA.elp',PATH));    
    EEG=pop_chanedit(EEG, 'lookup',sprintf('/Applications/eeglab13_6_5b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'));
    %re-reference to average
    EEG = pop_reref( EEG, []);

%%

 msec_matrix = dlmread(sprintf('%sData/OC/%s_OC.matrix', PATH, fn), '\t', 1, 1); % Namen einfügen
 
 EOG=EEG.data(61:64, :, :); % save EOG to put back after correction
 EEG.data = msec_matrix*EEG.data;
 EEG = pop_select(EEG,'channel',1:65);
 
 %% retrieve original EOG data for later use
 EEG.data(61:64,:,:) = EOG;
 
 
 ECG = EEG;
 %% Filtering 
 
 % We might not want to filter for time frequency analyses, however for ERP
 % stuff that's better!
 %data, lowcut, hicut
    %high-cut
    EEG = pop_eegfiltnew(EEG, [], 40); %0.5
    
    
    
 %%   save clean data for ERP analyses
  
    %save set
    pop_saveset(EEG, 'filename',sprintf('%sData/EEG_clean/%s_cleanNF',PATH,fn));

    
    %% retrieve original EEG
    
%     EEG=ECG;
%     
%     %filter low-cut only for TF analyses
%     EEG = pop_eegfiltnew(EEG, 0.5, []); 
%     pop_saveset(EEG, 'filename',sprintf('%sData/EEG_clean/%s_clean_TFA',PATH,fn));

end
%%
%chanlocs=EEG.chanlocs;
%save(sprintf('/Users/Romy/Dropbox (Brown)/ShenhavLab/EEG_ressources/Experiments/BASB_EEG/Data/Export/chanlocs.mat'), 'chanlocs'); 