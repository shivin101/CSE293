function params = TPM_Simulation_Parameters(filename,opt_type, condition)
%
% params = TPM_Simulation_Parameters(opt_type, condition)
%  
% This function sets the default parameters for a few sets of imaging
% conditions
% 
%   - filename          - Parameter output filename
%
%   - opt_type          - Optics type. Defaults to standard
%          'standard'   - Gaussian illumination of back aperture
%          'bessel'     - Bessel beam illumination (Lu et al 2017)
%          'stefo'      - Temporally focused beam (Prevedel et al 2016)
%          'vtwins'     - vTwINS illumination
%
%   - condition         - Experimental or prep condition
%          'standard'   - Default prep (GCaMP6f transgenic mouse cortex)
%          'gcamp7'     - GCaMP7f indicator replacing GCaMP6
%          'sparse'     - Only 10% of neurons are labeled
%          'nuclear'    - All neurons are only nuclearly labeled
%
% 2017 - Adam Charles and Alex Song

%% Process inputs

if(nargin<2)||isempty(opt_type)
  opt_type = 'standard';
end

if(nargin<3)||isempty(condition)
  condition = 'standard';
end

%% Setup parameters 
vol_params   = [];
vasc_params  = [];
neur_params  = [];
dend_params  = [];
axon_params  = [];
bg_params    = [];
psf_params   = [];
spike_opts   = [];
tpm_params   = [];
noise_params = [];

tpm_params.pavg         = 40;                                              % Power in units of mW

switch opt_type
  case 'lowNA'
    psf_params.NA       = 0.2;                                             % Numerical aperture of PSF
    psf_params.psf_sz   = [20 20 80];
  case 'standard'
    psf_params.NA       = 0.6;                                             % Numerical aperture of PSF
    psf_params.psf_sz   = [20 20 50];
  case 'bessel'
    psf_params = getDefaultPSFParams('bessel');    
    psf_params.psf_sz   = [20 20 80];
    tpm_params.pavg = 120;    
  case 'stefo'
    psf_params = getDefaultPSFParams('temporal-focusing');    
    psf_params.psf_sz   = [20 20 50];
    tpm_params.pavg = 400;
  case 'vtwins'
    psf_params = getDefaultPSFParams('vtwins');
    psf_params.psf_sz   = [20 20 80];
  otherwise
    error('Given optics type is invalid')
end

psf_params.objNA    = 0.8;                                                 % Numerical aperture of PSF
psf_params.zernikeWt = [0 0 0 0 0 0 0 0 0 0 0.06];                         % In units of wavelength, a small amount of spherical aberration and astigmatism added as uncorrected "system" aberrations

switch condition
  case 'standard'
  case 'gcamp7'
    spike_opts.prot = 'GCaMP7';
  case 'sparse'
    spike_opts.p_off = 0.95;
  case 'nuclear'
    scan_params.nuc_label = 1;
  otherwise
    error('Given optics type is invalid')
end

vol_params.vol_sz       = [500,500,100];                                   % Volume size to sample (in microns)
vol_params.vol_depth    = 100;
neur_params.fluor_dist  = [0 0];
spike_opts.nt           = 5000;
spike_opts.rate         = 10e-2;
spike_opts.burst_mean   = 10;                                              % Mean spike rate
scan_params.verbose     = 2;                                               % Set scanning verbosity level
scan_params.motion      = 1;                                               % Set motion option (true/false)
scan_params.scan_buff   = 10;                                              % Buffer set aside for motion

vol_params   = check_vol_params(vol_params);                               % Check volume parameters
vasc_params  = check_vasc_params(vasc_params);                             % Make default set of vasculature parameters
neur_params  = check_neur_params(neur_params);                             % Make default set of neuron parameters
dend_params  = check_dend_params(dend_params);                             % Make default set of dendrite parameters
axon_params  = check_axon_params(axon_params);                             % Make default set of axon parameters
bg_params    = check_bg_params(bg_params);                                 % Make default set of background parameters
psf_params   = check_psf_params(psf_params);                               % Check point spread function parameters
spike_opts   = check_spike_opts(spike_opts);                               % Check spike/fluorescence simulation parameters
tpm_params   = check_tpm_params(tpm_params);
scan_params  = check_scan_params(scan_params);
noise_params = check_noise_params(noise_params);                           % Make default noise parameter struct for missing elements

noise_params.sigscale = tpmSignalscale(tpm_params)*(spike_opts.dt)*(scan_params.sfrac^2)/(vol_params.vol_sz(1)*vol_params.vol_sz(2)*(vol_params.vres^2));


%% Output parameters
params.vol_params      = vol_params;
params.vasc_params     = vasc_params;
params.neur_params     = neur_params;
params.dend_params     = dend_params;
params.axon_params     = axon_params;
params.bg_params       = bg_params;
params.psf_params      = psf_params;
params.spike_opts      = spike_opts;
params.tpm_params      = tpm_params;
params.scan_params     = scan_params;
params.noise_params    = noise_params;

if((nargin>=1)&&~isempty(filename))
  save(filename,'-struct','params')
end
