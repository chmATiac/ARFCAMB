# ARFCAMB
CAMB modified version to include ARF as new observable

A set of routines have been modified to include angular redshift fluctuations (ARF) as a new cosmological observable. 

ARF may be invoked by stating " source_type='arf' " when calling redshift windows (temporaritly only the "SplinedSourceWindow" routine works properly, but "GaussianSourceWindow" will be modified in the future). The ARF will thus correspond to a given redshift window, just as it would be the case for " source_type='counts' " or " source_type='lensing' ". 
The impact of photometric redshift errors can only be included when using the "SplinedSourceWindow" routine. The RMS of the photo-z errors may be introduced by including the flag "sigma_Errz". The PDF of the photo-z errors is assumed to be Gaussian, but a (truncated) Lorentzian may be used if the flag "i_Lorentz=1" is added as an argument in the "SplinedSourceWindow" function. The Lorentzian is truncated at +/- 4 sigma_Errz, beyond which a Gaussian shape is adopted. Example in python: 


zS=np.arange(0.4,1.8,0.001) ; dndzS=(zS/0.6)**-0.5 

pars.SourceWindows = [ SplinedSourceWindow(bias=1.,source_type='arf',dlog10Ndm=0.0, z=zS, \
                        W=dndzS\*np.exp(-0.5*(zS-1.)\*\*2/0.01\*\*2), bias_z=np.sqrt(1+zS), sigma_Errz=0.01], i_Lorentz=0 )


If no explicit input is given for sigma_Errz and i_Lorentz, then the values of sigma_Errz=1e-5 and i_Lorentz=0 (Gaussian case) are taken as default values.

This modified version of CAMB also includes local "f_NL" non-Gaussianity parameter and the elliptical Non Linear Intrinsic Allignment model (eNLIA) for weak lensing. 
These are included within the standard set of CAMB parameters, and their defaul values are:

AIA = 1.72

epsIA = -0.41

betaIA = 2.17

for the eNL IA model, and 

f_NL = 0.0

for the local non-Gaussianity parameter.


