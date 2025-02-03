%% Second-level PEB analyses - Juliette Lanskey 2024
%  1. Compare baseline group differences (controls vs patients)
%  2. Compare longitudinal group differences (patients; baseline vs
%  follow-up)

%% Set up environment
clearvars
E = cmm_environment;

scr= E.scr;
anaL= E.anaL;
anaB= E.anaB;

load([anaB filesep 'Bsubs.mat']); 
Bbmr=load([anaB filesep 'RCM_BMC_BMA']); 
maxind=6;
BGCM=spm_dcm_load(Bbmr.RCM(:,maxind));

Lbmr=load([anaL filesep 'RCM_BMC_BMA.mat']);
LGCM=spm_dcm_load(Lbmr.RCM(:,maxind));

agecov = 1; %whether to include age as a PEB covariate in Baseline PEB
dayscov = 0; %whether to include days as a PEB covariate in longitudinal PEB


%% Define AMPA v GABA v NMDA fields
%--------------------------------------------------------------------------
field{1} = {'T(:,1)', 'A', 'B'};
field{2} = {'T(:,2)'};
field{3} = {'T(:,3)', 'AN', 'BN'};

labels = {'AMPA', 'GABA', 'NMDA'};

%% Baseline Group comparison
% specify PEB model
X=zeros(length(Bsubs),2);
X(:,1)=ones; %  mean

con = find(contains(Bsubs, 'C'));
pat = find(contains(Bsubs, 'P'));

X(con,2)=0;
X_labels={'mean','cons>pats'};
X(pat,2)=1;

if agecov==1
    
    load([scr 'BLage.mat']);
    for ss=1:length(Bsubs)
        idx=find(contains(BLage(:,1),Bsubs(ss)));
        Xage(ss,1) = BLage(idx,1);
        Xage(ss,2) = BLage(idx,2);
    end
    
    for ss=1:length(Bsubs)
        if strcmp(Xage{ss,1}, Bsubs{ss})
            X(ss,3)=Xage{ss,2};
        else
            disp('subject orders do not match')
        end
    end
    X_labels={'mean','cons>pats','age'};
    X(:,3:end)=zscore(X(:,3:end));
    
end

M=struct(); %a structure to specify the PEB settings
M.Q = 'all'; % between-subject variability will be estimated for each connection
M.X=X; %design matrix
M.Xnames = X_labels;

%% Run PEB recursively over the set of reduced models defined in 'C'
%--------------------------------------------------------------------------
for f = 1:length(field)
    PEB{f} = spm_dcm_peb(BGCM, M, field{f});
    F(f)   = PEB{f}.F;
end

%% Plot Bayesian model comparison of all PEB models
%--------------------------------------------------------------------------
plot_bmc_F_jl(F, labels);

if agecov == 1
    if ~exist([scr '/figures/PEB_comparison_B_age.png'], 'file')
        exportgraphics(gcf, [scr '/figures/PEB_comparison_B_age.png'], 'Resolution', '720');
    end
else 
    if ~exist([scr '/figures/PEB_comparison_B.png'], 'file')
        exportgraphics(gcf, [scr '/figures/PEB_comparison_B.png'], 'Resolution', '720');
    end
end

% look at effect on parameters for winning model
win=find(F==max(F));
f=win;
spm_dcm_peb_review_fig_jl(PEB{f}, 2, 0.95, 1);

%% Bayesian model comparison and averaging over symmetric space
RMA{f} = spm_dcm_peb_bmc(PEB{f});

Mtemp=spm_perm_mtx(10); Mtemp=full(Mtemp);
np = length(PEB{1}.Pnames);
Mall=zeros(size(Mtemp,1), np);
Mall(:, 1:8) = Mtemp(:, 1:8);
Mall(:, 11) = Mtemp(:, 9);
Mall(:, 15) = Mtemp(:, 10);

count=0;
for ss=1:size(Mall, 1)
    LHS = [Mall(ss,1:4), Mall(ss, 11)];
    RHS = [Mall(ss,5:8), Mall(ss, 15)];
    
    if all(LHS == RHS)
        count=count+1;
        sym_idx(count) = ss;
    end
end

Msymm=Mall(sym_idx,:);
RMAmodel = spm_dcm_peb_bmc(PEB{f},  Msymm);

spm_dcm_peb_review_fig_jl(RMAmodel, 2, 0.95, 2);
xticknames=xticklabels_jl_peb(PEB{f}, BGCM{1});
xticknames(3, 1:8)={'TN LA1', 'TN LSTG', 'TN LIFG', 'TN L IPC', 'TN RA1', 'TN RSTG', 'TN RIFG', 'TN RIPC'};
xticklabels([xticknames{3,:}])

if agecov == 1
    
    if ~exist([scr '/figures/group_NMDA_age.png'], 'file')
        exportgraphics(gcf, [scr '/figures/group_NMDA_age.png'], 'Resolution', '720');
    end
    
    if ~exist([scr '/figures/age_NMDA_age.png'], 'file')
        spm_dcm_peb_review_fig_jl(RMAmodel, 3, 0.95, 2);
        xticklabels([xticknames{3,:}])
        exportgraphics(gcf, [scr '/figures/age_NMDA_age.png'], 'Resolution', '720');
    end

else
    if ~exist([scr '/figures/group_NMDA.png'], 'file')
        exportgraphics(gcf, [scr '/figures/group_NMDA.png'], 'Resolution', '720');
    end
end

close all

%% Longitudinal PEB analysis
load([anaL '/Lsubs.mat']);
LX=zeros(2*length(Lsubs),2);
LX(:,1)=ones; %first collumn to ones to model overall mean

LX(1:length(Lsubs),2)=zeros;

if dayscov == 1
    load([scr '/days_af-bl.mat']); %should load into a variable called 'days'

    for ss=1:length(Lsubs)
        idx=find(contains(days(:,1),Lsubs(ss)));
        LX(length(Lsubs)+ss,3) = days{idx,2}/365.2425;
    end
    LX(length(Lsubs)+1:length(Lsubs)*2,3) = zscore(LX(length(Lsubs)+1:length(Lsubs)*2,3));
    LX(length(Lsubs)+1:end,2) = 1;
    LX_labels={'mean','BL->AF', 'days'};
    
else
    LX(length(Lsubs)+1:end,2) = 1;
    LX_labels={'mean','BL->AF'};
    
end

LM=struct(); %a structure to specify the PEB settings
LM.Q = 'all'; % between-subject variability will be estimated for each connection
LM.X=LX; %design matrix
LM.Xnames = LX_labels;

%% Comparing AMPA v GABA v NMDA for longitudinal analysis
%--------------------------------------------------------------------------
for f = 1:length(field)
    PEBL{f} = spm_dcm_peb(LGCM, LM, field{f});
    FL(f)   = PEBL{f}.F;
end

plot_bmc_F_jl(FL, labels);

if dayscov == 1
    if ~exist([scr '/figures/PEB_comparison_F_days.png'])
        exportgraphics(gcf, [scr '/figures/PEB_comparison_F_days.png'], 'Resolution', '720');
    end
else
    if ~exist([scr '/figures/PEB_comparison_F.png'])
        exportgraphics(gcf, [scr '/figures/PEB_comparison_F.png'], 'Resolution', '720');
    end
end

RMALm=spm_dcm_peb_bmc(PEBL{3}, Msymm);
spm_dcm_peb_review_fig_jl(RMALm, 2, 0.95, 2)
xticklabels([xticknames{3,:}])


if dayscov == 1
    if ~exist([scr '/figures/AF_NMDA_days.png'])
        exportgraphics(gcf, [scr '/figures/AF_NMDA_days.png'], 'Resolution', '720');
    end
    if ~exist([scr '/figures/days_NMDA_days.png'])
        spm_dcm_peb_review_fig_jl(RMALm, 3, 0.95, 2)
        xticklabels([xticknames{3,:}])
        exportgraphics(gcf, [scr '/figures/days_NMDA_days.png'], 'Resolution', '720');
    end
else
    if ~exist([scr '/figures/PEB_NMDA.png'])
        exportgraphics(gcf, [scr '/figures/PEB_NMDA.png'], 'Resolution', '720');
    end
end
