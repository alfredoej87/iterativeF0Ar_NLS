# iterativeF0Ar_NLS

%%Needs to previously download the fast F0 NLS estimator from https://github.com/jkjaer/fastF0Nls (i.e., the Maximum Likelihood Estimator under a white Gaussian noise assumption). 
%%The ICASSP 2020 results can be obtained from running the file expIterativeforICASSPresults.m (need to add the directory of fastF0Nls) 
%%The .pev and .pes files contain the audio samples of Keele database and the peak correlation values used for F0 estimation (assumed here as a "ground truth") 
%%Note that in babble noise conditions (and nonstationary noise, e.g., restaurant noise) a better improvement is achieved with using 256 noise spectral envelopes instead of the 16 for the ICASSP submission, when applying the initial pre-whitening (i.e., pre-processing) step. %%Even if some form of pre-processing (e.g, pre-whitening or speech enhancement) benefits non-parametric pitch estimators, a better accuracy from the introduced approach is observed, specially at low iSNRs, even if the resulting pitch estimates are not smoothed. 

%%Some results of other databases will follow as soon as possible. 
