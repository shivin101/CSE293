function psf = applyTemporalFocusing(psf,length,dz,offset)

% psf = applyTemporalFocusing(psf,length,dz,offset)
%
% This function applies temporal focusing restriction on the PSF by
% assuming the constraint on the axial profiles follows 1/sqrt(1+(dz/zR)^2)
% relationship, where dz is the offset axially and zR is the constrained
% width, and is related to the FWHM by a factor of 2*sqrt(3)
% This approximates the theoretical valeus from Durst 2006 (Optics Express)
% 
% Inputs:
%     psf     - 3D psf matrix (with axial scale dz)
%     length  - length of FWHM axially from temporal focusing (um)
%     dz      - step size in z per voxel
%     offset  - offset from center of psf matrix in z (in um)
% 
% Outputs:
%     psf     - temporally focused psf
%
% 2017 - Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin<4)
  offset = 0;
end
 
zR = length/(dz*2*sqrt(3));
z0 = ceil(size(psf,3)/2);

zprofile = 1./sqrt(1+((z0-(1:size(psf,3))+(offset*dz))/zR).^2);
psf = bsxfun(@times,reshape(zprofile,1,1,[]),psf);
