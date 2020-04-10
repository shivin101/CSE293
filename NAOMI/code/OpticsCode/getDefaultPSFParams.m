function psf_params = getDefaultPSFParams(psf_type)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch psf_type
case 'gaussian'
    psf_params.NA       = 0.5; 
case 'vtwins'
    psf_params.type     = 'vtwins';
    psf_params.psf_sz   = [20 80 80];
    psf_params.theta    = 30;                                              % angle between the two halves of the vTwINS PSF (degrees) (double value from midline)
    psf_params.sepdist  = 40;                                              % average separation distance between two halves of vTwINS PSF (um)
    psf_params.length   = 50;                                              % length of PSF (FWHM), used for NA;
    psf_params.propcrop = 0;                                               % crop edges during optical propagation flag (default 1)
    TMP = check_psf_params(psf_params);
    psf_params.n        = TMP.n;
    psf_params.lambda   = TMP.lambda;
    psf_params.beamNA = psf_params.n*sqrt(1-(0.626*psf_params.lambda/(psf_params.length*psf_params.n)-1)^2);
    psf_params.NA     = sin(atan(tan(asin(psf_params.beamNA))+tand(psf_params.theta/2)));
    psf_params.objNA  = psf_params.NA;
    case 'vtwins_bessel'
    psf_params.type       = 'vtwins';
    psf_params.vtwinstype = 'bessel';
    psf_params.psf_sz     = [20 80 80];
    psf_params.theta      = 30;                                            % angle between the two halves of the vTwINS PSF (degrees) (double value from midline)
    psf_params.sepdist    = 40;                                            % average separation distance between two halves of vTwINS PSF (um)
    psf_params.length     = 40;                                            % length of PSF (FWHM), used for NA;
    psf_params.propcrop   = 0;                                             % crop edges during optical propagation flag (default 1)
    psf_params.beamNA     = 0.3;
    psf_params.NA         = sin(atan(tan(asin(psf_params.beamNA))+tand(psf_params.theta/2)));
    psf_params.zernikeWt  = 0;
    psf_params.objNA      = psf_params.NA;
case 'bessel'
    psf_params.type      = 'bessel';
    psf_params.psf_sz    = [20 20 80];
    psf_params.NA        = 0.4;                                            % sets Bessel beam diameter (1.629*lambda*n/(2*pi*NA))
    psf_params.length    = 60;                                             % sets Bessel beam length
    psf_params.objNA     = 0.6;
    psf_params.zernikeWt = 0;
case 'temporal-focusing'
    psf_params.psf_sz   = [20 20 60];
    psf_params.length   = 10;                                              % 10um FWHM axial resolution
    psf_params.objNA    = 0.5;
    psf_params.NA       = 0.085;                                           % 5um FWHM lateral resolution
    psf_params.propcrop = 0;                                               % crop edges during optical propagation flag (default 1)
    psf_params.scaling  = 'temporal-focusing';
case 'three_photon'
    psf_params.scaling = 'three-photon';
end

psf_params   = check_psf_params(psf_params);                               % Check point spread function parameters
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
