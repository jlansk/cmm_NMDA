%% DCM using the cmm_NMDA model to model Alzheimer's disease progression
%  1. Inverts the full DCM for each individual and check fits
%  2. Runs a bayesian model reduction over model space

%% Set up environment
clearvars
E = cmm_environment;

scr= E.scr; 
anaL= E.anaL;
anaB= E.anaB;
raw= E.raw;

load([scr '/AFsubs.mat']);
subjects=AFsubs;
rawlast = '/AF/mmn/ffmraeMaffffdtsss.mat'; % file that has had all NTAD preprocessing steps applied)

%% Set up  model space
%create models from my model space with different A and C matrices
[DCMa] = dcm_cmm_gen_model_space_v1();

% specify rest of parameters for the full DCM model
for ss=1:length(subjects)
    dcms{ss} = dcm_cmm_gen_v1(DCMa{end}, raw, anaL, subjects{ss}, rawlast);
end

dcms=dcms';
GCM=spm_dcm_load(dcms);

%% Invert first-level DCMs
%  subject that do not converge go into the 'failed' variable

if invert == 1
    
    failed={};
    parfor ss=1:length(subjects)
        try
            spm_dcm_erp(GCM{ss});
        catch
            failed{ss}={strcat('failed_', num2str(ss))};
            continue;
        end
    end
    
    fail = find(~cellfun(@isempty,failed));
    failed_subs = subjects(fail);
    save([anaL '/failedsubs.mat'], 'failed_subs') 
end

%% Remove failed or bad fits
% From the variable 'failed_subs' created during the inversion and from
% seeing which have NANs in the H matrix (as shown through above plots) we
% remove those failed subjects
BLbmr = load([anaB filesep 'RCM_BMC_BMA.mat']);

if ~exist('maxind', 'var')
    maxind=6;
end

if ~exist('failed_subs', 'var')
    load([anaB '/failedsubs.mat'])
end


count=0;
for ss = 1:length(subjects)
    try
        DCM =spm_dcm_load([anaL filesep 'DCM_' subjects{ss} '_full.mat']);
        DCM = DCM{1};
        if any(isnan([DCM.H{2}(:); DCM.H{1}(:)]))
            count=count+1
            nan_subs{count}=subjects{ss}
        end
    end
end

subjects_orig=subjects;
clearvars subjects

count=0;
for ss=1:length(subjects_orig)
    if ~any(find(contains(subjects_orig{ss},[failed_subs, nan_subs])))
        count=count+1;
        
        subjects{count}=subjects_orig{ss};
    end
end

BGCM=spm_dcm_load(BLbmr.RCM(:,maxind));

for f=1:length(BGCM)
    Bfiles{f}=BGCM{f}.name;
end

Bsubs=extractBetween(Bfiles, 'DCM_', '_');
Lsubs = subjects(ismember(subjects, Bsubs)); %remove any subs not in BL who are in AF (failed to converge)
save([anaL '/Lsubs'], 'Lsubs');

%% Plot model fits (line graphs)
for ss = 1:length(subjects)
    
    DCM =spm_dcm_load([anaL filesep 'DCM_' subjects{ss} '_full.mat']);
    DCM = DCM{1};
    reps=1;
    
    subplot(4,15,ss)
    
    for c = 1
        plot(DCM.H{c}(:,1) + DCM.R{c}(:,1), 'color', [0.4 0.4 0.4], 'Linewidth', 1.5); hold on
        plot(DCM.H{c}(:,1), 'color', [0.3    0.8    0.6510], 'Linewidth', 1.5);
        
        obs1(ss,:)=DCM.H{c}(:,1) + DCM.R{c}(:,1);
        prd1(ss,:)=DCM.H{c}(:,1);
        
        cortemp=corrcoef(DCM.H{c}(:,1) + DCM.R{c}(:,1), DCM.H{c}(:,1));
        corB.s(ss, 1)=subjects(ss);
        corB.d(ss, 1)=cortemp(2);
        
        ylim([-7 4]); yticks([-7 4]);
        xlim([0 125]); xticks([25 75]);
        %xlabel(DCM.name(end-9:end-5))
       
        if strcmp(DCM.name(end-9), 'P')
            xlabel('Patient');
        elseif strcmp(DCM.name(end-9), 'C')
            xlabel('Control');
        end
        set(gcf, 'color', 'w');
        set(gcf, 'Position', [100 + 400*(reps-1) 100 400 800]);
    end
    box off
    
end

legend({'DEV(pred)', 'DEV(obs)'});

if ~exist([scr '/figures/Afits.png'])
    exportgraphics(gcf, [scr '/figures/Afits.png'], 'Resolution', 720) %only works in 2020+
end

%% get free energy of nested models for bayesian model comparison ...
[DCMa] = dcm_cmm_gen_model_space_v1();

% BASELINE: For each subject, generate reduced submodels to run BMR over
%--------------------------------------------------------------------------
ss = 1;

for s = 1:length(Lsubs)
    
    count = 1;
    
    
    for m = 1:(length(DCMa)-1)
        clear DCM
        DCM = dcm_cmm_gen_v1(DCMa{m}, raw, anaB, Lsubs{s}, rawlast);
        P{ss,count} = DCM;
        count = count + 1;
    end
    
    % Add the full model that has been inverted to the end
    %--------------------------------------------------------------------------
    FCM{s} = spm_dcm_load([anaB filesep 'DCM_' Lsubs{s} '_full.mat']);
    P{ss,count} = FCM{s}{1};
    
    ss = ss+1;
    
end

% FOLLOW UP
%--------------------------------------------------------------------------
%Generate reduced submodels from the full model for each subject

for s = 1:length(Lsubs)
    
    count = 1;
    
    %Add sub-models with exponential and phasic repetition effect
    %--------------------------------------------------------------------------
    for m = 1:(length(DCMa)-1)
        clear DCM
        DCM = dcm_cmm_05_gen(DCMa{m}, raw, anaL, Lsubs{s}, rawlast);
        P{ss,count} = DCM;
        count = count + 1;
    end
    
    
    % Add the full model that has been inverted to the end
    %--------------------------------------------------------------------------
    FCM{s} = spm_dcm_load([anaL filesep 'DCM_' Lsubs{s} '_full.mat']);
    P{ss,count} = FCM{s}{1};
    
    ss = ss+1;
end

%% Run BMR

[RCM, BMC, BMA] = spm_dcm_bmr(P);
save([anaL filesep 'RCM_BMC_BMA'], 'RCM', 'BMC', 'BMA');