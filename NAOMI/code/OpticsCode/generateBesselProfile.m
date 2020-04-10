function Uout = generateBesselProfile(X,Y,rad,width,k,fl,type,aper,offset)

% Uout = generateBesselProfile(X,Y,rad,aper,k,fl,offset)
%
% This function generates a Gaussian back aperture intensity profile with a
% fixed aperture size with a focused phase fl. The inputs are
%
% - X           - lateral position X [m]
% - Y           - lateral position X [m]
% - rad         - radius [m]
% - width       - width of annulus (FWHM) on back aperture [m]
% - k           - optical wavenumber [rad/m]
% - fl          - focal lenth of lens [m]
% - type        - one of 'gaussian', 'uniform', 'posax', 'negax' for shape
%                 of profile. Gaussian is symmetric, uniform is flat across
%                 the width, posax and negax are x*gauss(x) and are the
%                 profile formed by a positive or negative strength axicon
% - aper        - aperture distance [m]
% - offset      - offset of Gaussian position [m]
%
% The output is
% - Uout        - scalar field for a gaussian beam with apodization 1
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<9)
  offset = [0 0];
end
if(nargin<8)
  aper = inf;
end
if(nargin<7)
  type = 'gaussian';
end

X2 = X-offset(1);
Y2 = Y-offset(2);
rho = sqrt(X2.^2+Y2.^2);
if(strcmp(type,'gaussian'))
  width = width/(2*sqrt(log(2)));
  Uout = exp(-((rho-rad)/width).^2);  
elseif(strcmp(type,'uniform'))
  Uout = ((rho-rad)<(width/2)).*((rho-rad)>-(width/2));    
elseif(strcmp(type,'posax'))
  width = width/1.133;
  Uout = (rho-rad+width/sqrt(2)).*exp(-((rho-rad+width/sqrt(2))/width).^2).*(rho-rad+width/sqrt(2)>0);  
elseif(strcmp(type,'negax'))
  width = width/1.133;
  Uout = (rad-rho+width/sqrt(2)).*exp(-((rad-rho+width/sqrt(2))/width).^2).*(rad-rho+width/sqrt(2)>0);  
else
  warning('Specified type of Bessel beam unavailable, defaulting to Gaussian')  
end
  
rho2 = X.^2+Y.^2;
Uout = Uout.*(rho2<(aper^2)); % Fixed aperture
Uout = Uout.*exp(-1i*k/(2*fl)*rho2);                                       % Apply ideal phase
