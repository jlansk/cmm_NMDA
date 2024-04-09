%% Baseline - to be continued...
%maybe do controls in isolation and then pats combined. maybe leave this as
%is for now until I discuss with James
% Plot Bayesian model comparison results
%=========================================================

%% Set up environment
clearvars
E = cmm_environment;
anaB=E.anaB;
anaL=E.anaL;

Bbmr = load([anaB filesep 'RCM_BMC_BMA']);
BMC=Bbmr.BMC;

clearvars F Fm dF Fms Fsort dF maxind
figure
F = [BMC.F]';
Fma = mean(F,1);
Fm=Fma(1:6);
Fms = Fm - min(Fm);

Fsort = sort(Fm);
dF    = Fsort(end) - Fsort(end-1);
maxind= find(Fm==max(Fm)); nextind=find(Fm==Fsort(end-1));

subplot(2,1,1)
b = bar(Fms);
b(1).FaceColor= [0.4 0.4 0.4];
ylim([0 1.42*10^4])
title(['Relative model evidence for deviant effect across individuals: dF = ' num2str(dF)]);
ylabel('Free energy');% yticks([0 60000])
box off

Pm = softmax(Fm');
subplot(2,1,2)
b = bar(Pm)
b.FaceColor= [0.4 0.4 0.4];
box off
ylabel('Poseterior model probability'); yticks([0 1]);

set(gcf, 'color', 'white')
exportgraphics(gcf, [scr '/figures/IPC_BMC.png'], 'Resolution', '720');


%% Follow up


% 
% %% Plot results of Bayesian model comparison FREE ENERGY
% %%==========================================================================
% 
% 
if plot2 ==1
    
    
    %% Plot Bayesian model comparison results
    %%=========================================================  
    clearvars F Fm dF Fms Fsort dF maxind
    figure   
    AFbmr=load([anaL filesep 'RCM_BMC_BMA_BLAFonlysubs.mat']);
    F = [AFbmr.BMC.F]';

    load([anaL '/Lsubs.mat']);
   % if alltogether ==1
        clearvars F
        BFbmr=load([anaB filesep 'RCM_BMC_BMA_dec22.mat'])
        FB = [BFbmr.BMC.F]'
        FL = [AFbmr.BMC.F]'
        FLB=[FB; FL(length(Lsubs)+1:length(Lsubs)*2,:)];
        F = FLB;
 %   end
    
    Fma = mean(F,1);
    Fm=Fma(1:6)
    plot_bmc_F_jl(Fm, num2cell(1:1:6));
    
  %  exportgraphics(gcf, [anaL '/L_IPC_BMC.png'], 'Resolution', '720');

    
    % Controls alone
    load([anaB '/Bsubs.mat']);
    con = find(contains(Bsubs, 'C'));
    Fc = FB(con, :)
    Fmc = mean(Fc,1);
    Fmc=Fmc(1:6)
    plot_bmc_F_jl(Fmc, num2cell(1:1:6));
    title('controls only')
    
    % Patients baseline
    pat = find(contains(Bsubs, 'P'));
    Fp = FB(pat, :)
    Fmp = mean(Fp,1);
    Fmp=Fmp(1:6)
    plot_bmc_F_jl(Fmp, num2cell(1:1:6));
    title('patients baseline')
    
    % Patients follow-up
    pat = find(contains(Bsubs, 'P'));
    Fp = FL(length(Lsubs)+1:end, :)
    Fmp = mean(Fp,1);
    Fmp=Fmp(1:6)
    plot_bmc_F_jl(Fmp, num2cell(1:1:6));
    title('patients followup')
    
    % All patients baseline and follow-up
    pat = find(contains(Bsubs, 'P'));
    Fp = FLB(pat(1):end, :)
    Fmp = mean(Fp,1);
    Fmp=Fmp(1:6)
    plot_bmc_F_jl(Fmp, num2cell(1:1:6));
    title('patients all baseline, followup')
    
end