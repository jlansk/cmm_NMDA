function E = cmm_environment()

%% Set up environment
E.Fspm = '/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7771/';
addpath(E.Fspm);
spm('defaults', 'EEG');

E.scr='/imaging/rowe/users/jl01/meg/dcm_cmm/scripts/to_publish/';
addpath(genpath(E.scr));

E.raw='/imaging/rowe/users/jl01/meg/dcm_cmc/meg_data/'; % where the preprocessed data is kept

E.anaB =   '/imaging/rowe/users/jl01/meg/dcm_cmm/reps0-5/BL_Feb23/full/noA3'; % Baseline analysis folder
E.anaL =  '/imaging/rowe/users/jl01/meg/dcm_cmm/reps0-5/AF_Feb23/full/noA3'; % Longitudinal analysis folder

end
