function psf_params = check_psf_params(psf_params)

% psf_params = check_psf_params(psf_params)
%  
% This function checks the elements of the struct psf_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - psf_params        - Struct contaning the parameters for the PSF
%          .NA          - Numerical aperture of Gaussian beam
%          .n           - Refractive index of propagation material
%          .n_diff      - Shift in refractive index from vessels to tissue
%          .lambda      - Two-photon excitation wavelength (um)
%          .obj_fl      - Objective focal length (mm)
%          .ss          - Subsampling factor for fresnel propagation
%          .sampling    - Spatial sampling for tissue occlusion mask
%          .psf_sz      - Default two-photon PSF size simulated (um)
%          .prop_sz     - Fresnel propagation length outside of volume (um)
%          .blur        - PSF lateral blurring (um)
%          .scatter_sz  - Scattering object sizes (um), column vector
%          .scatter_wt  - Scattering object weights, column vector
%          .zernikeWt   - Microscope aberration weights (Zernike basis)
%
% 2017 - Adam Charles and Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(psf_params)
    clear psf_params
    psf_params = struct;
end

if (~isfield(psf_params,'NA'))||isempty(psf_params.NA)                     % Default excitation numerical aperture
    psf_params.NA = 0.6;
end
if (~isfield(psf_params,'objNA'))||isempty(psf_params.objNA)               % Default objective numerical aperture
    psf_params.objNA = 0.8;
end
if (~isfield(psf_params,'n'))||isempty(psf_params.n)                       % Default index of refraction in tissue
    psf_params.n = 1.35;
end
if (~isfield(psf_params,'n_diff'))||isempty(psf_params.n_diff)             % Default shift in index of refraction from vessels to tissue
    psf_params.n_diff = 0.02;
end
if (~isfield(psf_params,'lambda'))||isempty(psf_params.lambda)             % Default two-photon excitation wavelength (microns)
    psf_params.lambda = 0.92;
end
if (~isfield(psf_params,'obj_fl'))||isempty(psf_params.obj_fl)             % Default objective focal length (mm)
    psf_params.obj_fl = 4.5;
end
if (~isfield(psf_params,'ss'))||isempty(psf_params.ss)                     % Default subsampling factor for fresnel propagation (from volume voxel size)
    psf_params.ss = 2;
end
if (~isfield(psf_params,'sampling'))||isempty(psf_params.sampling)         % Default spatial sampling for tissue occlusion mask
    psf_params.sampling = 50;
end
if (~isfield(psf_params,'psf_sz'))||isempty(psf_params.psf_sz)             % Default two-photon PSF size simulated (microns)
    psf_params.psf_sz = [20 20 50];
end
if (~isfield(psf_params,'prop_sz'))||isempty(psf_params.prop_sz)           % Default fresnel propagation length outside of volume (microns)
    psf_params.prop_sz = 10;
end
if (~isfield(psf_params,'blur'))||isempty(psf_params.blur)                 % Default psf lateral blurring (microns)
    psf_params.blur = 3;
end
if (~isfield(psf_params,'scatter_sz'))                                     % Default scattering object sizes (microns), column vector
%     psf_params.scatter_sz = [0.21 1.00 3.41 13.11 92.01]';
      psf_params.scatter_sz = [0.51 1.56 4.52 14.78]';
end
if (~isfield(psf_params,'scatter_wt'))                                     % Default scattering object weights, column vector
%     psf_params.scatter_wt = 1.65*[1.0 0.11 0.020 0.0033 0.00016]';
    psf_params.scatter_wt = [0.57 0.29 0.19 0.15]';
end
if (~isfield(psf_params,'zernikeWt'))||isempty(psf_params.zernikeWt)       % Default microscope aberration weights (Zernike basis)
    psf_params.zernikeWt = [0 0 0 0 0.1 0 0 0 0 0 0.12];                    % In units of wavelength, a small amount of spherical aberration and astigmatism added as uncorrected "system" aberrations
end
if (~isfield(psf_params,'taillength'))||isempty(psf_params.taillength)     % Distance from edge of PSF_sz to estimate tailweight (um) 
    psf_params.taillength = 50;                                           
end
if (~isfield(psf_params,'type'))||isempty(psf_params.type)                 % Default PSF type ('gaussian','vtwins','bessel')
    psf_params.type = 'gaussian';                                    
end
if (~isfield(psf_params,'scaling'))||isempty(psf_params.scaling)           % Default PSF scaling type ('two-photon','three-photon','temporal-focusing')
    psf_params.scaling = 'two-photon';                                    
end
if (~isfield(psf_params,'hemoabs'))||isempty(psf_params.hemoabs)           % Hemoglobin absorbance scaling factor
    psf_params.hemoabs = 0.00674*log(10);                                  % Default assumes 150mg/ml Hb, 64500 g/mol Hb, 2.9 (abs/um)/(mol/L) in units of abs/um. Absorbance calculated from Scott Prahl's Hb curve and eGFP emission spectrum
end
if (~isfield(psf_params,'propcrop'))||isempty(psf_params.propcrop)         % Flag to crop scanned beam during optical propagation (default true)
    psf_params.propcrop = true;                                            
end
if (~isfield(psf_params,'fastmask'))||isempty(psf_params.fastmask)         % Flag to crop scanned beam during optical propagation (default true)
    psf_params.fastmask    = true;
    psf_params.FM.sampling = 10;
    psf_params.FM.fineSamp = 2;
    psf_params.FM.ss       = 1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%