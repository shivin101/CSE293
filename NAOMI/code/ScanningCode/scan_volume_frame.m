function [clean_img, z_off] = scan_volume_frame(scan_vol, neur_act, scan_params, z_off)

% mov = scan_volume(neur_vol, PSF, neur_act, scan_params, varargin)
%
% Scan a volume and create a simulated two-photon calcium imaging movie.
% The inputs to this function are: 
%   neur_vol     - A struct that contains the output of the volume
%                  generating code. 
%   PSF          - A struct representing the point-spread function of the
%                  microscopy setup
%   neur_act     - Struct containing the simulated fluorescence activity 
%                  for each of the K neurons for all nt time-steps
%    .soma       - Kx(nt) matrix where each row is the activity of each 
%                  cell at the soma 
%    .dend       - Kx(nt) matrix where each row is the activity of each 
%                  cell at the dendrites [can be optional]
%    .bg         - N_bgx(nt) matrix where each row is the activity of each 
%                  background/neuropil component [can be optional]
%   scan_params  - Struct contaning the scanning parameters. Includes
%    .scan_avg   - Sampling rate of the scanning in terms of how many
%                  granular pixels to scan into one pixel (default = 2)
%    .motion     - True/false option for whether or not to simulate
%                  motion while performing scanning (default = true)
%    .scan_buff  - Number of granular pixels to keep as a buffer from the
%                  edge of the volume (default = 10)
%    .verbose    - Level of verbosity in the output during the volume
%                  generation. Can be 0,1,2. 0 = no text updates, 1 = 
%                  some text outputs. 2 = detailed text outputs (default
%                  = 1)
%   noise_params - Struct containing parameters for the noise model
%    .mu         - Mean measurement increase per photon (default = 0.5)
%    .mu0        - Electronics offset (default = 0)
%    .sigma      - Variance increase per photon (default = 0.3)
%    .sigma0     - Electronics base noise variance (default = 0.1)
%
% The outputs of this function is: 
%   mov         - A 3D array where each 2D slice is one frame in the video
%                 recroding of the simulated microscopy data
%
% 2016 - Adam Charles and Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup parameters

% scan_buff = scan_params.scan_buff;                                         % Get scanning side buffer number (how much to stay away from the edges of the volume)
% mot_opt   = scan_params.motion;                                            % Get scanning motion option (true or false) to decide whether to simulate tissue-FOV motion
scan_avg  = scan_params.scan_avg;                                          % Get scanning stepping amount (how many sub-resolution steps between each sample point in the FOV)
sfrac     = scan_params.sfrac;                                             % Subsampling factor

N1     = scan_params.vol_sz(1);                                          % Get length dimension
N2     = scan_params.vol_sz(2);                                          % Get width dimension
N3     = scan_params.vol_sz(3);                                          % Get depth dimension

Np1    = scan_params.psf_sz(1);
Np2    = scan_params.psf_sz(2);
Np3    = scan_params.psf_sz(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize movie and temporary volume
% x_loc = 1;
% y_loc = 1;
if(nargin<4)
  z_off = 0;
else
  zmaxdiff = 2;
  z_off = z_off+(rand<0.005)-(rand<0.005);
  z_off = max(-zmaxdiff,z_off);
  z_off = min(z_off,zmaxdiff);
end

z_loc = floor(0.5*(N3-Np3))+z_off;


% if (~isfield(scan_params,'zoffset'))||(isempty(scan_params.zoffset))
%   z_base = floor(0.5*(N3-Np3));
% else
%   z_base = floor(0.5*(N3-Np3))+scan_params.zoffset;
% end
% z_loc  = z_base;                                                           % Set scan depth
% x_loc  = floor(scan_buff+1);                                               % Set initial scan location in x
% y_loc  = floor(scan_buff+1);                                               % Set initial scan location in y

% if mot_opt                                                                 % If motion is requested, initialize the motion simulation variables
%     zmaxdiff = 2;%5;
%     d_stps   = [-1,1,zeros(1,5)];                                          % Set the vector of steps to sample from when choosing how much to move the volume relative to the FOV at each frame
%     d_stpsZ  = [-1,1,zeros(1,100)];                                        % Set the vector of steps to sample from when choosing how much to move the volume (axially) relative to the FOV at each frame
%     d_stps2  = -3:3;                                                       % Set the list of admissible steps (in granular resolution steps) that the volume can move relative to the FOV due to larger motion jumps
%     p_jump   = 0.05;                                                       % Set the probability of a bigger jump in the FOV
%     maxshear = 1/200;                                                      % fraction of max shearing that can be from fast y motion
% else                                                                       % ----- 
%     zmaxdiff = 0;
%     d_stps   = [0,0,0];                                                    %   | 
%     d_stpsZ  = [0,0,0];                                                    % If no motion is requested, set all the motion simulation variables to zero (this will stop all motion)
%     d_stps2  = [0,0,0];                                                    %   |
%     p_jump   = 0;                                                          %   |
%     maxshear = 0;                                                          % -----
% end

%% Setup scanning volume
if(~isstruct(neur_act))
  act = neur_act;
  neur_act.soma = act;
  neur_act.dend = act;
  neur_act.bg = act;
  clear act
end
cutoff = 1e-2;

if (isfield(scan_params,'nuc_label'))&&(scan_params.nuc_label>=1);
    nuc_act = single(neur_act.soma);
    nuc_act(nuc_act<cutoff) = 0;
    
    neur_act.nuc = neur_act.soma;
    neur_act.soma = 0*neur_act.soma;
    neur_act.dend = 0*neur_act.dend;
    neur_act.bg = 0*neur_act.bg;
end

soma_act = single(neur_act.soma);
soma_act(soma_act<cutoff) = 0;

dend_act = single(neur_act.dend);
dend_act(dend_act<cutoff) = 0;

bg_act = single(neur_act.bg);
bg_act(bg_act<cutoff) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iteratively scan the volume

% if rand(1) > p_jump                                                % Randomly choose whether to jump
%   x_loc = min(max(1,x_loc+randsample(d_stps2,1)), 2*scan_buff+1);% Update random jump for x-axis (slow) scan start
%   y_loc = min(max(1,y_loc+randsample(d_stps2,1)), 2*scan_buff+1);% Update random jump for y-axis (fast) scan start
% end
% 
% x_pos = min(max(1,x_loc+randsample(d_stps, 1)), 2*scan_buff+1);    % Update random drift for x-axis (slow) scan start
% y_pos = min(max(1,y_loc+randsample(d_stps, 1)), 2*scan_buff+1);    % Update random drift for y-axis (fast) scan start
% z_loc = min(max(z_base-zmaxdiff,z_loc+...
%   randsample(d_stpsZ, 1)), z_base+zmaxdiff);          % Update random drift for depth scan start
% z_loc = min(max(1,z_loc+randsample(d_stps, 1)), N3-Np3+1);         % Update random drift for depth scan start
% 
% y_shr = [zeros(randsample(floor(2*N1/5),1),1); ...
%   linspace(0,1,round(rand*3*N1/5))'*(2*(rand-0.5))*maxshear*N1];    % Create shearing vector
% y_shr = cat(1,y_shr, y_shr(end)*ones(max(0,N1-numel(y_shr)),1));   % Make sure vector is the full size (N1)
% y_off = vec(min(max(1,y_pos+y_shr + ...
%   randsample(d_stps, N1, true)'),2*scan_buff+1));    % Randomly select offset for each row (+/-0.5um), in addition to the shearing vector

TMPvol = zeros(N1,N2,N3,'single');                                 % Reset the temporary volume - It's an order of magnitde faster to scrap and start over

for ll = 1:size(soma_act,1)                                        % For large volumes, mex functions decrease runtimes a lot
  if(soma_act(ll)>0)&&(~isempty(scan_vol.somaVol{ll,2}))
    array_SubSubTest(TMPvol,scan_vol.somaVol{ll,1},...
      scan_vol.somaVol{ll,2},soma_act(ll));
  end
end

if(isfield(scan_params,'nuc_label'))&&(scan_params.nuc_label>=1);
  for ll = 1:size(neur_vol.gp_nuc,1)
    if(nuc_act(ll)>0)&&(~isempty(scan_vol.nucVol{ll,2}))
      array_SubSubTest(TMPvol,scan_vol.nucVol{ll,1},...
        scan_vol.nucVol{ll,2},nuc_act(ll));
    end
  end
end

for ll = 1:size(dend_act,1)                                        % For large volumes, mex functions decrease runtimes a lot
  if (dend_act(ll)>0)&&(~isempty(scan_vol.dendVol{ll,1}))
    array_SubSubTest(TMPvol,scan_vol.dendVol{ll,1},...
      scan_vol.dendVol{ll,2},dend_act(ll));
  end
end

for ll = 1:size(bg_act,1)
  if(~isempty(scan_vol.axonVol{ll,1}) && bg_act(ll)>0)
    array_SubModTest(TMPvol, scan_vol.axonVol{ll,1}, ...
      scan_vol.axonVol{ll,2},bg_act(ll)); 
  end
end

clean_img = (1/(2*sfrac^2))*single_scan(TMPvol(:,:,...
  z_loc:(z_loc+Np3-1)),[Np1 Np2 Np3], scan_vol.freq_psf, scan_avg, true); % Scan a single frame, normalization for pixel size and scanning code

if(~isempty(scan_vol.g_blur))                                                % Keep clean_img with no additive blur
  clean_img = clean_img + conv2(clean_img,scan_vol.g_blur,'same');        % Add in additional blurring for dimmer areas
end

if(~isempty(scan_vol.psfT))
  if(isempty(1:z_loc-1))
    top_img = blurredBackComp2(TMPvol,1:size(TMPvol,3),...
      scan_vol.psfT.freq_psf,scan_vol.psfT.weight, scan_vol.top_mask, 1, []);
  else
    top_img = blurredBackComp2(TMPvol,1:z_loc-1,scan_vol.psfT.freq_psf,...
      scan_vol.psfT.weight, scan_vol.top_mask, 1, scan_vol.psfT.psfZ);
  end
  if(isempty(z_loc+Np3:size(TMPvol,3)))
    bot_img = blurredBackComp2(TMPvol,1:size(TMPvol,3),...
      scan_vol.psfB.freq_psf,scan_vol.psfB.weight, scan_vol.bot_mask, 1, []);
  else
    bot_img = blurredBackComp2(TMPvol,z_loc+Np3:size(TMPvol,3),...
      scan_vol.psfB.freq_psf,scan_vol.psfB.weight, scan_vol.bot_mask, 1, scan_vol.psfB.psfZ);
  end
  top_img = top_img/(sfrac^2);
  bot_img = bot_img/(sfrac^2);
  
  clean_img = clean_img + top_img + bot_img;
end

% clean_img   = imgSubRowShift(clean_img, scan_buff, x_pos, ...
%   round(y_off)); % Sub-select a portion of an image with row shifting

clean_img   = conv2(clean_img, ones(sfrac,sfrac), 'same');         % Sum up nearby pixels
clean_img   = clean_img(1:sfrac:end,1:sfrac:end);                  % Subsample for actual sampling rate

% motion = [x_loc;y_loc;z_loc];                          % If needed store the current scan location to the output array

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
