function scan_ideal(neur_vol, PSF_struct, neur_act, scan_params, varargin)

% mov = scan_ideal(neur_vol, PSF, neur_act, scan_params, varargin)
%
% Scan a volume and create ideal scanned components.
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
%% Parse inputs

if ~isfield(PSF_struct, 'psf')
    error('Must provide PSF to scan!')
else
    PSF    = PSF_struct.psf;                                               % Extract point spread function from PSF struct
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

if ~isfield(PSF_struct, 'g_blur')
    g_blur = [];                                                            % No additional blurring if transversal blur function not supplied
else
    g_blur = PSF_struct.blur;                                              % Extract point transversal blur function from PSF struct
end

if nargin > 4
    noise_params = varargin{1};                                            % Check if noise model parameters are provided
else
    noise_params = struct;                                                 % If not, initialize a struct to be filled in
end

scan_params  = check_scan_params(scan_params);                             % Check the scanning parameter struct for missing elements
noise_params = check_noise_params(noise_params);                           % Check the noise parameter struct for missing elements

if ~isstruct(neur_act)                                                     % Make sure that neur_act is a struct
    if size(neur_act,3)==1
        neur_act = struct('soma',neur_act);                                % Can be a Kx(nt)x1 array if only soma activity is provided
    elseif size(neur_act,3)==2
        neur_act = struct('soma',neur_act(:,:,1),'dend',neur_act(:,:,2));  % Can be a Kx(nt)x2 array if soma and dendrite activity is provided
    end
end
if (~isfield(neur_act,'dend'))||(isempty(neur_act.dend))
    neur_act.dend = neur_act.soma;                                         % Default to dendrites having the same activity as somas, if no dendrite activity is provided
end
if (~isfield(neur_act,'bg'))||(isempty(neur_act.bg))
    neur_act.bg = zeros(1,size(neur_act,2),'single');                      % Default to no background if background is not provided
end
if (isfield(scan_params,'nuc_label'))&&(scan_params.nuc_label>=1);
    neur_act.nuc = neur_act.soma;
    neur_act.soma = 0*neur_act.soma;
    neur_act.dend = 0*neur_act.dend;
    neur_act.bg = 0*neur_act.bg;
end

if (~isfield(scan_params,'fsimCleanPath'))
  scan_params.fsimCleanPath = [];
end

if (~isfield(scan_params,'saveBlocksize'))||(isempty(scan_params.saveBlocksize))
  scan_params.saveBlocksize = 1000;
end

if (~isfield(scan_params,'movout'))||(isempty(scan_params.movout))
  scan_params.movout = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize data

n0 = min(neur_act.soma,[],2);                                              % Get minimum value for each neuron's soma activity
n1 = min(neur_act.dend,[],2);                                              % Get minimum value for each neuron's dendrite activity
n2 = min(neur_act.bg,[],2);                                                % Get minimum value for each background component's activity

gIdxs = (n0+n1+n2)>0;

neur_base.soma = zeros(size(neur_act.soma,1),sum(gIdxs),'single');         % Initialize the somatic component profiles
neur_base.dend = zeros(size(neur_act.dend,1),sum(gIdxs),'single');         % Initialize the dendritic component profiles
neur_base.bg   = zeros(size(neur_act.bg,1),sum(gIdxs),'single');           % Initialize the background component profiles

neur_base.soma(sub2ind(size(neur_base.soma),find(gIdxs)',1:sum(gIdxs))) = n0(gIdxs);
neur_base.dend(sub2ind(size(neur_base.dend),find(gIdxs)',1:sum(gIdxs))) = n1(gIdxs);
neur_base.bg(  sub2ind(size(neur_base.bg),  find(gIdxs)',1:sum(gIdxs))) = n2(gIdxs);

neur_act = neur_base;
clear neur_base
if(scan_params.scan_avg>1)
    warning('Warning: Scan averages greater than 1 not supported for more than 1 offset')
end
scan_avg = scan_params.scan_avg;                                           % Extract scanning average value
sfrac    = scan_params.sfrac;                                              % Ectract subsampling factor

if (~isfield(scan_params,'vol_sz'))||(isempty(scan_params.vol_sz))
  N1 = size(neur_vol.neur_vol,1);                                          % Get length dimension
  N2 = size(neur_vol.neur_vol,2);                                          % Get width dimension
  N3 = size(neur_vol.neur_vol,3);                                          % Get depth dimension
else
  N1 = scan_params.vol_sz(1);                                              % Get length dimension
  N2 = scan_params.vol_sz(2);                                              % Get width dimension
  N3 = scan_params.vol_sz(3);                                              % Get depth dimension
end
Nt = size(neur_act.soma,2);                                                % Get number of time steps

if ~isstruct(PSF)                                                          % If the PSF is not a struct (i.e., should be an array)
    [Np1,Np2,Np3] = size(PSF);                                             % Get size of the PSF 
elseif isfield(PSF,'left')
    [Np1,Np2,Np3] = size(PSF.left);                                        % Reserved for vTwINS testing: Probably best not to play with this unless you know what you're doing
else
    error('Unknown input configuration!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform some error checking

if (N1 < Np1)||(N2 < Np2)
    error('PSF extent is bigger than the volume!')                         % Check that the PSF fits inside the volume (transversally)
end
if (N3 < Np3)
    error('PSF depth is larger than the volume depth!')                    % Check that the PSF fits inside the volume (axially)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize movie and temporary volume

TMPvol = zeros(N1,N2,N3,'single');                                         % Initialize temporary volume
if (~isfield(scan_params,'zoffset'))||(isempty(scan_params.zoffset))
  z_base = floor(0.5*(N3-Np3));
else
  z_base = floor(0.5*(N3-Np3))+scan_params.zoffset;
end

%%%%%%%%%%%%%%% TO DO: MAKE INTO A FIELD
zmaxdiff = 2; 
%%%%%%%%%%%%%%%

if (isfield(PSF_struct,'psfT'))&&(~isempty(PSF_struct.psfT))
    psfT = PSF_struct.psfT;
    psfB = PSF_struct.psfB;
    psfT.mask = psfT.mask/mean(psfT.mask(:));
    psfB.mask = psfB.mask/mean(psfB.mask(:));
    psfT.convmask = psfT.convmask/sum(psfT.convmask(:));
    psfB.convmask = psfB.convmask/sum(psfB.convmask(:));
    psfT.freq_psf = psf_fft(size(TMPvol), psfT.convmask);  
    psfB.freq_psf = psf_fft(size(TMPvol), psfB.convmask);  
else
    psfT = [];
    psfB = [];
end

freq_psf = psf_fft(size(TMPvol), PSF, scan_avg);                           % Pre-calculate the FFT of the PSF for faster scanning
sigscale = noise_params.sigscale;                                          % scaling factor for signal
somaVol  = cell(size(neur_vol.gp_vals,1),2);                               % --- 
dendVol  = cell(size(neur_vol.gp_vals,1),2);                               %  |
for i = 1:size(neur_vol.gp_vals,1)                                         %  |
    somaVol{i,1} = neur_vol.gp_vals{i,1}(neur_vol.gp_vals{i,3});           % Pre-alocate separate soma/dendrite indexing for faster activity modulation in the volume
    dendVol{i,1} = neur_vol.gp_vals{i,1}(~neur_vol.gp_vals{i,3});          %  |
    if(isempty(t_mask))                                                    %  |
        somaVol{i,2} = neur_vol.gp_vals{i,2}(neur_vol.gp_vals{i,3});       %  |
        dendVol{i,2} = neur_vol.gp_vals{i,2}(~neur_vol.gp_vals{i,3});      %  |
    else                                                                   %  |
        somaVol{i,2} = neur_vol.gp_vals{i,2}(neur_vol.gp_vals{i,3}).* ...  %  |
                                    t_mask(max(1,mod(somaVol{i,1},N1*N2)));%  |
        dendVol{i,2} = neur_vol.gp_vals{i,2}(~neur_vol.gp_vals{i,3}).* ... %  |
                                    t_mask(max(1,mod(dendVol{i,1},N1*N2)));%  |      
    end                                                                    %  |
end                                                                        % ---
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup scanning volume

cutoff = 1e-2;

soma_act = single(neur_act.soma);
soma_min = min(soma_act,[],2);
soma_act = bsxfun(@minus,soma_act,soma_min);
soma_act(soma_act<cutoff) = 0;

dend_act = single(neur_act.dend);
dend_min = min(dend_act,[],2);
dend_act = bsxfun(@minus,dend_act,dend_min);
dend_act(dend_act<cutoff) = 0;

bg_act = single(neur_act.bg);
bg_min = min(bg_act,[],2);
bg_act = bsxfun(@minus,bg_act,bg_min);
bg_act(bg_act<cutoff) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iteratively scan the volume

if (~isempty(scan_params.fsimCleanPath))
    tifLinkFsimClean = cell(2*zmaxdiff+1,1);
    tagFsimClean = cell(2*zmaxdiff+1,1);
    fsimCleanPath = cell(2*zmaxdiff+1,1);
    for jj = 1:(2*zmaxdiff+1)
        [path,name,~] = fileparts(scan_params.fsimCleanPath);
        fsimCleanPath{jj} = fullfile(path,sprintf([name '_%02d.tif'],jj));
        [tifLinkFsimClean{jj},tagFsimClean{jj}] = tifinitialize(fsimCleanPath{jj}, ...
              floor([(N1)/sfrac, (N2)/sfrac,]));
    end
end

if scan_params.verbose == 1
    fprintf('Scanning...')
elseif scan_params.verbose > 1
    fprintf('Scanning...\n')
end
if ~isstruct(PSF)
    for kk = 1:Nt                                                          % Iterate over time steps
        if scan_params.verbose >1
            tic
        end
        clear TMPvol; 
        TMPvol = zeros(N1,N2,N3,'single');                                 % Reset the temporary volume - It's an order of magnitde faster to scrap and start over                                                
        for ll = 1:size(soma_act,1)                                        % For large volumes, mex functions decrease runtimes a lot
            if(soma_act(ll,kk)>0)&&(~isempty(somaVol{ll,2}))
                array_SubSubTest(TMPvol,somaVol{ll,1},...
                                        somaVol{ll,2},soma_act(ll,kk));    % Iteratively add in each neuron's soma activity
            end
        end 
        if(isfield(scan_params,'nuc_label'))&&(scan_params.nuc_label>=1);
            for ll = 1:size(neur_vol.gp_nuc,1)
                if(nuc_act(ll,kk)>0)&&(~isempty(nucVol{ll,2}))
                    array_SubSubTest(TMPvol,nucVol{ll,1},...
                                        nucVol{ll,2},nuc_act(ll,kk));    % Iteratively add in each neuron's soma activity
                end            
            end
        end
        for ll = 1:size(dend_act,1)                                        % For large volumes, mex functions decrease runtimes a lot
            if (dend_act(ll,kk)>0)&&(~isempty(dendVol{ll,1}))
                array_SubSubTest(TMPvol,dendVol{ll,1},...
                                           dendVol{ll,2},dend_act(ll,kk)); % Iteratively add in each neuron's soma activity
            end
        end
        for ll = 1:size(bg_act,1)
          if(~isempty(axonVol{ll,1}) && bg_act(ll,kk)>0)
            array_SubModTest(TMPvol, axonVol{ll,1}, ...
                                           axonVol{ll,2},bg_act(ll,kk)); % Iteratively add in each background component's activity
          end
        end
        
        z_loc = floor(0.5*(N3-Np3));
        clean_img_all = (sigscale/(2*sfrac^2))*single_scan_stack(TMPvol,size(PSF), freq_psf, scan_avg, true, z_loc:(z_loc+Np3-1),-zmaxdiff:zmaxdiff);
          
        for jj = 1:(2*zmaxdiff+1)
          z_loc = jj-zmaxdiff-1+floor(0.5*(N3-Np3));          
          clean_img = clean_img_all(:,:,jj);
          if(~isempty(g_blur))                                                % Keep clean_img with no additive blur
            clean_img = clean_img + conv2(clean_img,g_blur,'same');        % Add in additional blurring for dimmer areas
          end
          if(~isempty(psfT))
            if(~isempty(t_mask))
                top_mask = 1./t_mask;
                bot_mask = 1./t_mask;
            else
                top_mask = ones(size(clean_img,1),size(clean_img,2),'single');
                bot_mask = ones(size(clean_img,1),size(clean_img,2),'single');
            end
            if(isfield(psfT,'mask'))
                top_mask = top_mask.*psfT.mask;
                bot_mask = bot_mask.*psfB.mask;
            end
            if(isempty(1:z_loc-1))
                top_img = blurredBackComp2(TMPvol,1:size(TMPvol,3),...
                                  psfT.freq_psf,psfT.weight, top_mask, 1, []);
            else
                top_img = blurredBackComp2(TMPvol,1:z_loc-1,psfT.freq_psf,...
                                         psfT.weight, top_mask, 1, psfT.psfZ);
            end
            if(isempty(z_loc+Np3:size(TMPvol,3)))
                bot_img = blurredBackComp2(TMPvol,1:size(TMPvol,3),...
                                  psfB.freq_psf,psfB.weight, bot_mask, 1, []);
            else
                bot_img = blurredBackComp2(TMPvol,z_loc+Np3:size(TMPvol,3),...
                           psfB.freq_psf,psfB.weight, bot_mask, 1, psfB.psfZ);
            end
            top_img = top_img*(sigscale/(sfrac^2));
            bot_img = bot_img*(sigscale/(sfrac^2));
            
            clean_img = clean_img + top_img + bot_img;
          end
          clean_img   = conv2(clean_img, ones(sfrac,sfrac), 'same');         % Sum up nearby pixels
          clean_img   = clean_img(1:sfrac:end,1:sfrac:end);                  % Subsample for actual sampling rate
          
          if (~isempty(scan_params.fsimCleanPath))
              tifLinkFsimClean{jj} = tifappend(tifLinkFsimClean{jj},clean_img, ...
                tagFsimClean{jj},scan_params.saveBlocksize,kk,fsimCleanPath{jj});
          end
        end
        if scan_params.verbose == 1
            fprintf('.');
        elseif scan_params.verbose >1
            Ttmp = toc;
            fprintf('Scanned frame %d (%f s)\n',kk,Ttmp)
        end
    end
elseif isfield(PSF,'left')                                                 % Check if separate beam v_twins configuration is being used
    error('Supported PSF configuration!')    
else
    error('Unknown input configuration!')
end

if (~isempty(scan_params.fsimCleanPath))
  for jj = 1:(2*zmaxdiff+1)
    tifLinkFsimClean{jj}.close();
  end
end


if scan_params.verbose >= 1
    fprintf('done.\n')
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
