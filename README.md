# cmm_NMDA
NMDA-channel kinetics of people living with Alzheimer’s disease
This code accompanies Lanskey et al. (2024).

The code runs on Matlab 2020b and uses SPM12 v7771, which is freely available.


## Analysis steps
The analysis is run in the following order:
1. Preprocessing, using the folder '1_preprocessing' from * [dcm_cmc_ntad](https://github.com/jlansk/dcm_cmc_ntad/tree/main/scripts/1_preprocessing)
2. Sensor space analysis, folder '2_sensor space'
3. First-level DCM inversion, folder '3_dcm'
4. Second-level Parameteric Empirical Bayes, folder 4_peb
