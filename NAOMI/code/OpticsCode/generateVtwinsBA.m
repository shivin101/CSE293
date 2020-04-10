function Uout2 = generateVtwinsBA(vol_params,psf_params)

% Uout2 = generateVtwinsBA(vol_params,psf_params)
%
% This function generates a Vtwins back aperture intensity profile with 
% the specified aberrations applied
%
%
% - vol_params          - Struct with parameters for the volume generation
%          .vol_sz      - 3-element vector with the size (in um) of the 
%                         volume to generate (default = 100x100x30um)
%          .vres        - resolution to simulate volume at (default = 2
%                         samples/um)
% - psf_params          - Struct contaning the parameters for the PSF
%          .NA          - Numerical aperture of Gaussian beam
%          .n           - Refractive index of propagation material
%          .lambda      - Two-photon excitation wavelength (um)
%          .obj_fl      - Objective focal length (mm)
%          .ss          - Subsampling factor for fresnel propagation
%          .sampling    - Spatial sampling for tissue occlusion mask
%          .zernikeDst  - Microscope aberration weights (Zernike basis) as
%                         a function of position
%
% - Uout2               - Output scalar field
%
% 2017 - Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fl = single(psf_params.obj_fl/1000);                                       % focal length [m]
D2 = single(1e-6*(1/vol_params.vres)/psf_params.ss);                       % observation grid spacing [m]
N = single(1e-6*(vol_params.vasc_sz(1:2)-vol_params.vol_sz(1:2))/D2);      % gridsize
D1 = single(max(gaussianBeamSize(psf_params,fl*1e6)/1e6)/min(N));          % source grid spacing [m]
nre = single(psf_params.n);                                                % immersion numerical index
rad = single(tan(asin(psf_params.objNA/nre))*fl);                          % source radius [m]
k = 2*nre*pi/single(psf_params.lambda*1e-6);                               % optical wavenumber [rad/m]
[X,Y] = meshgrid((-N(1)/2:N(1)/2-1)*D1,(-N(2)/2:N(2)/2-1)*D1);             

beamRad = single(tan(asin(psf_params.beamNA/nre))*fl);                     % half beam radius [m]
beamOffset = single(tand(psf_params.theta/2))*fl;                          % source radius [m]

if(~isfield(psf_params,'vtwinstype')||isempty(psf_params.vtwinstype))
  UoutL = generateGaussianProfile(X,Y,beamRad,rad,k,fl,[beamOffset 0]);
  UoutR = generateGaussianProfile(X,Y,beamRad,rad,k,fl,[-beamOffset 0]);
elseif(strcmp(psf_params.vtwinstype,'bessel'))
  width = 1.4e-10*tan(asin(psf_params.NA/nre))/(fl*1e-6*psf_params.length);
  if(~isfield(psf_params,'bprof')||isempty(psf_params.bprof))
    UoutL = generateBesselProfile(X,Y,beamRad,width,k,fl,'gaussian',rad,[beamOffset 0]);
    UoutR = generateBesselProfile(X,Y,beamRad,width,k,fl,'gaussian',rad,[-beamOffset 0]);
  else
    UoutL = generateBesselProfile(X,Y,beamRad,width,k,fl,psf_params.bprof,inf,[beamOffset 0]);
    UoutR = generateBesselProfile(X,Y,beamRad,width,k,fl,psf_params.bprof,inf,[-beamOffset 0]);
  end
end

tilt = [0 1 0]*psf_params.lambda*1e-6*psf_params.sepdist/8;
UoutL = applyZernike(UoutL,X/rad,Y/rad,k,tilt);
UoutR = applyZernike(UoutR,X/rad,Y/rad,k,-tilt);

imax = round(vol_params.vol_sz(1)/psf_params.sampling)+1;
jmax = round(vol_params.vol_sz(2)/psf_params.sampling)+1;

if(imax*jmax==1||(~isfield(psf_params,'zernikeDst'))||isempty(psf_params.zernikeDst))
  abb = generateZernike(psf_params);
  Uout2L = applyZernike(UoutL,X/rad,Y/rad,k,abb);
  Uout2R = applyZernike(UoutR,X/rad,Y/rad,k,abb);
else
  Uout2L = cell(imax,jmax);
  Uout2R = cell(imax,jmax);
  for i = 1:imax
    for j = 1:jmax
      abb = generateZernike(psf_params);
      Uout2L{i,j} = applyZernike(Uout2L,X,Y,k,abb);
      Uout2R{i,j} = applyZernike(Uout2R,X,Y,k,abb);
    end
  end
end
Uout2.left = Uout2L;
Uout2.right = Uout2R;