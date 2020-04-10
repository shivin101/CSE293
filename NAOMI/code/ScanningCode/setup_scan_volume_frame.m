function scan_vol = setup_scan_volume_frame(neur_vol,PSF_struct,scan_params)

scan_params  = check_scan_params(scan_params);                             % Check the scanning parameter struct for missing elements

scan_avg  = scan_params.scan_avg;                                          % Get scanning stepping amount (how many sub-resolution steps between each sample point in the FOV)
sfrac     = scan_params.sfrac;                                             % Subsampling factor

N1     = scan_params.vol_sz(1);                                          % Get length dimension
N2     = scan_params.vol_sz(2);                                          % Get width dimension
N3     = scan_params.vol_sz(3);                                          % Get depth dimension

if ~isfield(PSF_struct, 'psf')
    error('Must provide PSF to scan!')
else
    PSF    = PSF_struct.psf;                                               % Extract point spread function from PSF struct
end

if ~isstruct(PSF)                                                          % Get size of the PSF
    [Np1,Np2,Np3] = size(PSF);
else
    error('Unknown input configuration!')
end

if (N1 < Np1)||(N2 < Np2)
    error('PSF extent is bigger than the volume!')                         % Check that the PSF fits inside the volume (transversally)
end
if (N3 < Np3)
    error('PSF depth is larger than the volume depth!')                    % Check that the PSF fits inside the volume (axially)
end


if (isfield(PSF_struct,'psfT'))&&(~isempty(PSF_struct.psfT))
  psfT = PSF_struct.psfT;
  psfB = PSF_struct.psfB;
  psfT.mask = psfT.mask/mean(psfT.mask(:));
  psfB.mask = psfB.mask/mean(psfB.mask(:));
  psfT.convmask = psfT.convmask/sum(psfT.convmask(:));
  psfB.convmask = psfB.convmask/sum(psfB.convmask(:));
  psfT.freq_psf = psf_fft([N1 N2 N3], psfT.convmask);  
  psfB.freq_psf = psf_fft([N1 N2 N3], psfB.convmask);  
else
  psfT = [];
  psfB = [];
end

if ~isfield(PSF_struct, 'mask')
    t_mask = [];                                                           % No masking if mask is not supplied
else
    t_mask = PSF_struct.mask;                                              % Extract mask from PSF struct
    t_thresh = 1e-5;
    t_mask(t_mask<t_thresh) = t_thresh;
end
if isfield(PSF_struct, 'colmask')
    t_mask = t_mask.*PSF_struct.colmask;                                   % No masking if mask is not supplied
end

if(~isempty(t_mask))
  top_mask = 1./t_mask;
  bot_mask = 1./t_mask;
else
  top_mask = ones(N1/sfrac,N2/sfrac,'single');
  bot_mask = ones(N1/sfrac,N2/sfrac,'single');
end
if(isfield(psfT,'mask'))
  top_mask = top_mask.*psfT.mask;
  bot_mask = bot_mask.*psfB.mask;
end

scan_vol.top_mask = top_mask;
scan_vol.bot_mask = bot_mask;

clear top_mask bot_mask

scan_vol.psfT = psfT;
scan_vol.psfB = psfB;

clear psfT psfB

scan_vol.freq_psf = psf_fft([N1 N2 N3], PSF, scan_avg);                           % Pre-calculate the FFT of the PSF for faster scanning


if ~isfield(PSF_struct, 'g_blur')
    scan_vol.g_blur = [];                                                            % No additional blurring if transversal blur function not supplied
else
    scan_vol.g_blur = PSF_struct.blur;                                              % Extract point transversal blur function from PSF struct
end


somaVol  = cell(size(neur_vol.gp_vals,1),2);                               % --- 
dendVol  = cell(size(neur_vol.gp_vals,1),2);                               %  |
for i = 1:size(neur_vol.gp_vals,1)                                         %  |
    somaVol{i,1} = neur_vol.gp_vals{i,1}(neur_vol.gp_vals{i,3});           % Pre-alocate separate soma/dendrite indexing for faster activity modulation in the volume
    dendVol{i,1} = neur_vol.gp_vals{i,1}(~neur_vol.gp_vals{i,3});          %  |
    if(isempty(t_mask))
        somaVol{i,2} = neur_vol.gp_vals{i,2}(neur_vol.gp_vals{i,3});       %  |
        dendVol{i,2} = neur_vol.gp_vals{i,2}(~neur_vol.gp_vals{i,3});      %  |
    else
        somaVol{i,2} = neur_vol.gp_vals{i,2}(neur_vol.gp_vals{i,3}).* ...
                                    t_mask(max(1,mod(somaVol{i,1},N1*N2)));%  |
        dendVol{i,2} = neur_vol.gp_vals{i,2}(~neur_vol.gp_vals{i,3}).* ...
                                    t_mask(max(1,mod(dendVol{i,1},N1*N2)));%  |      
    end
end                                                                        % ---

scan_vol.somaVol = somaVol;
clear somaVol;
scan_vol.dendVol = dendVol;
clear dendVol;

if(~isempty(neur_vol.bg_proc))                                             % t_mask multiplied by gp_vals{i,2}
    if(isempty(t_mask))
      axonVol = neur_vol.bg_proc;
    else
      axonVol = cell(size(neur_vol.bg_proc));
      for i = 1:size(neur_vol.bg_proc,1)
        axonVol{i,1} = neur_vol.bg_proc{i,1};
        axonVol{i,2} = neur_vol.bg_proc{i,2}.* ...
             t_mask(max(1,mod(axonVol{i,1},N1*N2)));
      end
    end
else
    axonVol = [];
end

scan_vol.axonVol = axonVol;
clear axonVol;

if(isfield(scan_params,'nuc_label'))&&(scan_params.nuc_label>=1);
    if(isempty(t_mask))
      nucVol = neur_vol.gp_nuc;
    else
      nucVol = cell(size(neur_vol.gp_nuc));
      for i = 1:size(neur_vol.gp_nuc,1)
        nucVol{i,1} = neur_vol.gp_nuc{i,1};
        nucVol{i,2} = t_mask(max(1,mod(nucVol{i,1},N1*N2)));
      end
    end
  scan_vol.nucVol = nucVol;
  clear nucVol;
end

