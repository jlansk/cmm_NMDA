%% Sensor level analysis - Juliette Lanskey 2024
%This script runs the sensor space t-tests for the baseline
%(controls vs AD/MCI) and longitudinal (AD/MCI baseline vs follow up)
%comparisons and generates mismatch negativity waveform plots.

% set up variables
clearvars
E = cmm_environment;
raw = E.raw; %directory for MEG files
scr = E.scr; % directory of scripts
ana_dir = raw; % where the pre-processed MEG files are

load([scr filesep 'BLsubs.mat']); % 1x59 cell with baseline IDs
load([scr filesep 'AFsubs.mat']); % 1x33 cell with follow-up IDs

task='mmn';
session1 = 'BL';
session2 = 'AF';
lfile= 'bCPffmraeMaffffdtsss.mat'; % file with NTAD preprocessing pipeline and additional sensor-space processing (combined gradiometers, then baseline correction)

%% load baseline MMN amplitude
[mmn_pa_BL, mmnBL] = mmn_amp(BLsubs, ana_dir, 'BL', task, lfile, 0);

%% load follow-up MMN amplitude
[mmn_pa_AF, mmnAF] =mmn_amp(AFsubs, ana_dir, 'AF', task, lfile, 0);

%% statistics

%Baseline Group comparison
cons = find(contains(BLsubs, 'C'));
pats = find(contains(BLsubs, 'P'));

[h,p,ci,stats]=ttest2(mmn_pa_BL(cons,5),mmn_pa_BL(pats,5), 'tail','left');
cnanmean=nanmean(mmn_pa_BL(cons,1)); pnanmean=nanmean(mmn_pa_BL(pats,1));
BLstatistics(1,1)=stats.tstat; BLstatistics(2,1)=p; BLstatistics(3,1)=cnanmean; BLstatistics(4,1)=pnanmean;

%Longitudinal comparison
pats = find(contains(BLsubs, 'P'));

Lsubs = AFsubs(ismember(AFsubs, BLsubs));
LBsubs = find(contains(BLsubs, Lsubs));
LAsubs = find(contains(AFsubs, Lsubs));

[h,p,ci,stats]=ttest(mmn_pa_BL(LBsubs,5),mmn_pa_AF(LAsubs,5), 'tail','left');
bnanmean=nanmean(mmn_pa_BL(LBsubs,1)); fnanmean=nanmean(mmn_pa_AF(LAsubs,1));

AFstatistics(1,1)=stats.tstat; AFstatistics(2,1)=p; AFstatistics(3,1)=bnanmean; AFstatistics(4,1)=fnanmean;

%% plots
% MMN waveforms
BL_av = squeeze(nanmean(mmnBL(LBsubs,:,:),1))
BL_se= squeeze(nanstd(mmnBL(LBsubs,:,:),1))./sqrt(size(LBsubs,2));

AF_av = squeeze(nanmean(mmnAF(LAsubs,:,:),1))
AF_se = squeeze(nanstd(mmnAF(LAsubs,:,:),1))./sqrt(size(LAsubs,2))


con_av = squeeze(nanmean(mmnBL(cons,:,:),1))
con_se=squeeze(nanstd(mmnBL(cons,:,:),1))./sqrt(size(cons,2))%standard error

pat_av = squeeze(nanmean(mmnBL(pats,:,:),1))
pat_se = squeeze(nanstd(mmnBL(pats,:,:),1))./sqrt(size(pats,2))

% To save or not to save?

titlestr={'MMN_r5r0'};
choice = menu('Save figures?', 'Yes', 'No');


%% BL AD VS HC
close all

bh(1)=plot(con_av(:,5), 'Color', [1 0.42 0.16]); hold on; bh(2)=plot(pat_av(:,5), 'Color', [0 0.32 1]); hold on;
set(bh(1), 'linewidth', 3); set(bh(2), 'linewidth', 3);
boundedline([1:size(con_av)],con_av(:,5),con_se(:,5), 'transparency', 0.6, 'cmap', [1 0.55 0.16]); boundedline([1:size(pat_av)],pat_av(:,5),pat_se(:,5), 'transparency', 0.4, 'cmap', [0 0.1 0.75]); %see uisetcolor for colour options
xlim([0 250]); xticks([0 50 100 150 200 250]); xticklabels({'-100','0', '100', '200', '300', '400'}); xlabel('Time (ms)');
xlabel('Time (ms)'); ylabel('Mismatch response (fT/m)');
begin2=140; beg2=(begin2+100)/2;
finish2=160 ; fin2=(finish2+100)/2;
ylim([-1 0.1])
% patch([beg2,beg2, fin2, fin2, beg2], [0.1, -1.0, -1.0, 0.1, 0.1], 'black', 'EdgeColor', 'none', 'FaceColor', 'black', 'FaceAlpha', 0.1);
box off
set(gcf, 'color', 'w')

legend('controls', 'patients'); legend('Location', 'best')

if choice ==2 | choice ==0
    return;
elseif choice ==1
        exportgraphics(gcf, [scr filesep '/figures/CP_' titlestr{1} '.png'], 'Resolution', 720) %only works in 2020
end

%% longitudinal plot
close all
lh(1)=plot(BL_av(:, 5), 'Color', [0 0 0.6]); hold on; lh(2)=plot(AF_av(:, 5), 'Color', [0.05 0.6 0.65]); hold on;
set(lh(1), 'linewidth', 3); set(lh(2), 'linewidth', 3);
boundedline([1:size(BL_av)],BL_av(:,5),BL_se(:,5), 'transparency', 0.333, 'cmap', [0 0.1 0.7]); 
boundedline([1:size(AF_av)],AF_av(:,5),AF_se(:,5), 'transparency', 0.333, 'cmap',  [0.05 0.6 0.75]); %[0.05 0.8 0.85]see uisetcolor for colour options
xlim([0 250]); xticks([0 50 100 150 200 250]); xticklabels({'-100','0', '100', '200', '300', '400'}); xlabel('Time (ms)');
xlabel('Time (ms)'); ylabel('Mismatch response (fT/m)');
%  patch([beg2,beg2, fin2, fin2, beg2], [0.1, -1.0, -1.0, 0.1, 0.1], 'black', 'EdgeColor', 'none', 'FaceColor', 'black', 'FaceAlpha', 0.1);
begin2=140; beg2=(begin2+100)/2; 
finish2=160 ; fin2=(finish2+100)/2;
ylim([-1 0.1]); box off; set(gcf, 'color', 'w');
legend('baseline', 'follow-up'); legend('Location', 'best')

if choice ==2 | choice ==0
    return;
elseif choice ==1
        exportgraphics(gcf, [scr filesep '/figures/BF_' titlestr{1} '.png'], 'Resolution', 720) %only works in 2020 onwards  
end
    
