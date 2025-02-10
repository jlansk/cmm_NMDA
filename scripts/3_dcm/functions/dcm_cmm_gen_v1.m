function DCM = dcm_cmm_gen_v1(DCMa, raw, ana, sub, rawlast)

% Generating the DCM parameters
%==========================================================================

% DCMa is your model specification of A B and C matrices. 
% raw is the start of your file address. 
% sub is subject name
% rawlast is end of file address of data file to be inverted
% ana is your analysis directory

        
%--------------------------------------------------------------------------
%% ---- Set A and C matrices of DCM from DCMa       
DCM = DCMa; %specifies DCM.A, DCM.C

% Data filename
%--------------------------------------------------------------------------

DCM.xY.Dfile =  [raw filesep sub rawlast]; 


%% ----DCM model parameters ----
% Parameters and options used for setting up model
%--------------------------------------------------------------------------
DCM.options.analysis = 'ERP'; % analyze evoked responses
DCM.options.model    = 'CMM_NMDA'; %'CMM_NMDA'; % CMC, TFM
DCM.options.Tdcm(1)  = 0;     % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = 300;   % end of peri-stimulus time to be modelled
DCM.options.Nmodes   = 8;     % nr of modes for data selection
DCM.options.h        = 1;     % nr of DCT components
DCM.options.onset    = 60;    % selection of onset (prior mean)
DCM.options.dur     = 16;
DCM.options.D        = 1;     % downsampling
DCM.options.trials   = [1 6]; % index of ERPs within file

DCM.xY.modality='MEGPLANAR';
DCM.options.spatial  = 'ECD'; %LFP 

DCM.options.location = 1;     % optimising source location
DCM.options.han      = 1;     % applying hanning window
DCM.options.symmetry=1;

% Location priors for dipoles
%--------------------------------------------------------------------------
%coords from Lappe2013 for IPC
DCM.Lpos = [[-42; -22; 7] [ -61; -32; 8] [ -46; 20; 8] [-58; -27; 30] [ 46; -14; 8] [ 59; -25; 8] [ 46; 20; 8] [59; -41; 30]]; 
DCM.Sname  = {'left A1'; 'left STG'; 'left IFG'; 'left IPC'; 'right A1'; 'right STG'; 'right IFG'; 'right IPC'};
Nareas    = size(DCM.Lpos,2);

% Specify B matrix
%--------------------------------------------------------------------------
DCM.B{1} = DCM.A{1} + DCM.A{2} + DCM.A{3};  % model specification
DCM.B{1}=logical(DCM.B{1});

        
% Between trial effects
%--------------------------------------------------------------------------
% Trials defined as deviant standard
%% ---- Between-trial effects ----

DCM.xU.X(:,1) = [1; 0];
DCM.xU.name{1} = ('d-r5');
DCM.xU.dt = 0.002; 
        
% Define priors
%--------------------------------------------------------------------------
[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
DCM.M.pE = pE;
DCM.M.pC = pC;

DCM.name = [ana filesep 'DCM_' sub '_' DCMa.name];
       
end

