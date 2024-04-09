%% DCM using the cmm_NMDA model to test baseline differences (patients vs controls)
%  1. Inverts the full DCM for each individual and checks fits
%  2. Runs a bayesian model reduction to compare models to identify where
%  parietal node is connected to other regions


%% Set up environment
clearvars
E = cmm_environment;

addpath(genpath('/imaging/rowe/users/jl01/meg/dcm_cmc/scripts/to_publish/'))

%% Set up variables
scr = E.scr;

load([scr '/BLsubs.mat']);
subjects = BLsubs;
anaB = E.anaB; %character array to DCM files location
raw = E.raw;  % character aray to base file directory (e.g. raw='/juliette/cmm/meg_data')
rawlast = '/BL/mmn/ffmraeMaffffdtsss.mat'; % file that has had all NTAD preprocessing steps applied


%% Set up  model space

%create models from my model space with different A and C matrices
[DCMa] = dcm_cmm_gen_model_space_v1();


% specify rest of parameters for the full DCM model
for ss=1:length(subjects)
    dcms{ss} = dcm_cmm_gen_v1(DCMa{end}, raw, anaB, subjects{ss}, rawlast);
end

dcms=dcms';
GCM=spm_dcm_load(dcms);

%% Invert first-level DCMs
%  subject that do not converge go into the 'failed' variable

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
save([anaB '/failedsubs.mat'],'failed_subs')


%% Remove failed or bad fits
% From the variable 'failed_subs' created during the inversion and from
% seeing which have NANs in the H matrix we remove those failed subjects

if ~exist('failed_subs', 'var')
    load([anaB '/failedsubs.mat'])
end

count=0; nan_subs={};
for ss = 1:length(subjects)
    try
        DCM =spm_dcm_load([anaB filesep 'DCM_' subjects{ss} '_full.mat']);
        DCM = DCM{1};
        if any(isnan([DCM.H{2}(:); DCM.H{1}(:)]))
            count=count+1
            nan_subs{count}=subjects{ss}
        end
    end
end

count=0;
for ss=1:length(subjects)
    if ~any(find(contains(subjects{ss},[failed_subs, nan_subs])))
        count=count+1;
        Bsubs{count}=subjects{ss};
    end
end
save([anaB '/Bsubs.mat'], 'Bsubs');

%% Plot model fits
for ss = 1:length(Bsubs)
    
    DCM =spm_dcm_load([anaB filesep 'DCM_' Bsubs{ss} '_full.mat']);
    DCM = DCM{1};
    
    subplot(6,15,ss); reps=1;
    
    for c = 1
        plot(DCM.H{c}(:,1) + DCM.R{c}(:,1), 'color', [0.4 0.4 0.4], 'Linewidth', 1.5); hold on
        plot(DCM.H{c}(:,1), 'color', [0.3    0.8    0.6510], 'Linewidth', 1.5); hold on
        
        obs1(ss,:)=DCM.H{c}(:,1) + DCM.R{c}(:,1);
        prd1(ss,:)=DCM.H{c}(:,1);

        cortemp=corrcoef(DCM.H{c}(:,1) + DCM.R{c}(:,1), DCM.H{c}(:,1));
        corB.s(ss, 1)=subjects(ss);
        corB.d(ss, 1)=cortemp(2);
             
        ylim([-5 5]); yticks([-5 5]);
        xlim([0 125]); 
        xticks([25 75]);
        xticklabels({'50', '150'})
        %xlabel(extractBetween(DCM.name, 'DCM_', '_'))
        
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

if ~exist([scr '/figures/Bfits.png'])
    exportgraphics(gcf, [scr '/figures/Bfits.png'], 'Resolution', 720) %only works in 2020
end

%% Bayesian model reduction
if ~exist('Bsubs', 'var')
    load([anaB '/Bsubs.mat']) %
end

[DCMa] = dcm_cmm_gen_model_space_v1();

%   Generate reduced submodels from the full model for each subject
for ss = 1:length(Bsubs)
    count = 1;
    
    for m = 1:(length(DCMa)-1)
        clear DCM
        DCM = dcm_cmm_gen_v1(DCMa{m}, raw, anaB, Bsubs{ss}, rawlast);
        P{ss,count} = DCM;
        count = count + 1;
    end
    
    %    Add the full model that has been inverted to the end
    %   --------------------------------------------------------------------------
    FCM{ss} = spm_dcm_load([anaB filesep 'DCM_' Bsubs{ss} '_full.mat']);
    P{ss,count} = FCM{ss}{1};
end


%  Run BMR
[RCM, BMC, BMA] = spm_dcm_bmr(P);

save([anaB filesep 'RCM_BMC_BMA'], 'RCM', 'BMC', 'BMA');