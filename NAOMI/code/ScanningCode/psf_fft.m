function freq_psf = psf_fft(vol_sz, psf, varargin)

% scan_img = scan_volume(neur_vol, psf, varargin)
%
% Scan a 3D volume with a given point-spread function. This function takes
% in a neural volume (as generated by "simulate_neural_vol_v*.m"), a
% point-spread function, and a number of other optional scanning parameters
% in order to simulate the images resulting from a two-photon scan of the
% volume with the given PSF. Inputs to this function are:
%   - neur_vol - 3D volume where each voxel contains the fluorescence
%   - psf      - 3D array containing the intensity of the point-spread
%                function
%   - z_sub    - OPTIONAL speed-up parameter that scans multiple slices
%                simultaneously, reducing the number of convolutions needed
%                by a factor of 'z_sub'
%
% The output is
%   - scan_img - The scanned image - i.e. the return fluorescence from the
%                point-spread function with no noise (photon or electronic)
%
% 2016 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse Inputs

if nargin > 2
    z_sub = varargin{1};
else
    z_sub = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate and sum convolutions

if z_sub > 1                                                               % For time considerations, can pre-sum every few slices
    N_slce    = ceil(size(psf,3)/z_sub);                                     % Figure out how many slices after pre-summing
    psf2      = psf(:,:,1:z_sub:z_sub*N_slce);                             % Initialize the new point-spread function
    for kk = 2:z_sub
        slcs                = kk:z_sub:min(z_sub*N_slce, size(psf,3));       % Slices to add
        Nz                  = numel(slcs);                                 % Figure out how many slices in the kk^th sub-sampling (in case 
        psf2(:,:,1:Nz)      = psf2(:,:,1:Nz) + psf(:,:,slcs);              % Iteratively add in the PSF slices
    end
    sz = vol_sz + size(psf2) - 1;                                          % Get the sizes of the post-convolution array
    sz = nearest_small_prime(sz,7);                                        % Make sure the largest factor is no larger than 7 (for fast FFTs)
    freq_psf = fft(fft(psf2, sz(1), 1), sz(2), 2);                         % Get fft of psf
else
    try
        sz = vol_sz + size(psf) - 1;                                       % Get the sizes of the post-convolution array
    catch
        sz = vol_sz + [size(psf) 1] -1;
    end
    sz = nearest_small_prime(sz,7);                                        % Make sure the largest factor is no larger than 7 (for fast FFTs)
    freq_psf = fft(fft(psf, sz(1), 1), sz(2), 2);                          % Get fft of psf
end
 


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
