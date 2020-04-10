%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATE TPM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Two-photon simulation code. This script shows how to use the tools in
% this code package to simulate a neural volume, neural activity, and
% scanning via a simulated two-photon imaging setup. The resulting videos
% mimic the statistics and activity seen in TPM data and can be used in
% evaluating optical parameter set-ups as well as assessing analysis
% algorithms. This code function in five main modules: neural simulation,
% volume creation, activity simulation, optical set-up, and scanning. The
% main function(s) for each of these are:
%
%    neuron creation     - simulate_neural_volume
%    volume creation     - simulate_neural_volume
%    activity simulation - generateTimeTraces
%    optical set-up      - simulate_optical_propagation
%    scanning simulation - scan_volume
%
% Each function has instructions available using the "help" command, in
% addition to the use example in this file. Additional help is available by
% contacting the creators, Adam Charles (adamsc@princeton.edu and Alex
% Song. 
% 
% 2019 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Parameters
profile on                         
save_loc = 'E:\NAOMIvols\example_vol.mat';
vol_params.vol_sz   = [250,250,100];                                       % Volume size to sample (in microns)
vol_params.vol_depth = 100;

vol_params   = check_vol_params(vol_params);                               % Check volume parameters
vasc_params  = check_vasc_params([]);                                      % Make default set of vasculature parameters
neur_params  = check_neur_params([]);                                      % Make default set of neuron parameters
dend_params  = check_dend_params([]);                                      % Make default set of dendrite parameters
axon_params  = check_axon_params([]);                                      % Make default set of axon parameters
bg_params    = check_bg_params([]);                                        % Make default set of background parameters
noise_params = check_noise_params([]);                                     % Make default noise parameter struct for missing elements
psf_params   = check_psf_params([]);                                       % Check point spread function parameters
tpm_params   = check_tpm_params([]);

debug_opt            = 0;
vol_params.verbose   = 0;

debug_opt = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate volume
[~,vol_params,neur_params,vasc_params,dend_params,bg_params,axon_params] ...
    = simulate_neural_volume_lowram(vol_params, neur_params, vasc_params, ...
     dend_params, bg_params, axon_params, psf_params, save_loc, debug_opt, 0); % Draw a random volume - this takes the longest amound of time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reload volume

load(save_loc,'vol_params', 'neur_params','vasc_params', 'dend_params', 'bg_params', 'axon_params', 'psf_params');
vol_out = load(save_loc,'neur_ves_all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create PSF/optical simulation

psf_params.NA = 0.6;
PSF_struct = simulate_optical_propagation(vol_params,psf_params,vol_out);  % Create the point-spread function and mask for scanning
PSF_struct.mask = PSF_struct.mask.^1.5;
save(save_loc,'-append','PSF_struct','psf_params');
clear PSF_struct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reload the neural locations only

load(save_loc,'vol_params', 'neur_params','vasc_params', 'dend_params', 'bg_params', 'axon_params', 'psf_params');
vol_out = load(save_loc,'neur_locs');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate neural activity

clear spike_opts
spike_opts.nt   = 1000;
spike_opts.rate = 0.16;
spike_opts.K    = vol_params.N_neur+vol_params.N_den+vol_params.N_den2;    % Read off the number of neurons
spike_opts.scalevec = [0.8 1.3 1.6];
spike_opts      = check_spike_opts(spike_opts);
act_rate = spike_opts.dt;
spike_opts.dt   = 1/100;

clear cal_params
cal_params.sat_type  = 'Ca_DE';                                     % Set the calcium dynamics to 'single' mode (simpler dynamics)
cal_params.prot_type = spike_opts.prot;                                           % Extract the type of protein to simulate
cal_params = check_cal_params(cal_params, cal_params.prot_type);                      % Check param struct and fill in with defaults

mod_vals = single(expression_variation(spike_opts.K, spike_opts.p_off, spike_opts.min_mod));    % Create a set of modulatory factors for each cell
resampStruct = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up scanning simulation

clear scan_params
scan_params.verbose   = 2;                                                 % Set scanning verbosity level
scan_params.sfrac     = 1;
scan_params.motion    = 1;                                                 % Set motion option (true/false)
scan_params.scan_buff = 5;
scan_params           = check_scan_params(scan_params);

tpm_params.pavg       = 40;                                                % power in units of mW

vol_out = load(save_loc,'bg_proc','gp_vals','gp_nuc');
load(save_loc,'PSF_struct')

scan_params.vol_sz = vol_params.vol_sz*vol_params.vres;
scan_params.psf_sz = size(PSF_struct.psf);

scan_vol = setup_scan_volume_frame(vol_out,PSF_struct,scan_params);
clear vol_out PSF_struct
load(save_loc,'neur_locs');
%%
[path,name,~] = fileparts(save_loc);
datapath = fullfile(path,[name '_data.h5']);
h5create(datapath,'/spikes',[spike_opts.K ceil((spike_opts.nt+1)*act_rate/spike_opts.dt)]);
h5create(datapath,'/act',[spike_opts.K spike_opts.nt]);
h5create(datapath,'/motion',[3 spike_opts.nt]);
h5create(datapath,'/cleanimg',[scan_params.vol_sz(1:2)/scan_params.sfrac spike_opts.nt]);
nt = spike_opts.nt*(act_rate/spike_opts.dt);
k = 0;
z_off = 0;
for i = 1:spike_opts.nt
  countprint(i);
  numpts = floor(i*act_rate/spike_opts.dt)-floor((i-1)*act_rate/spike_opts.dt);
  spikes = zeros(spike_opts.K,numpts,'single');
  act = zeros(spike_opts.K,numpts,'single');
  for j = 1:numpts
    [resampStruct,spikes(:,j),act(:,j)]  = generateNextTimePoint(spike_opts,cal_params,neur_locs,mod_vals,resampStruct); % Generate time traces using AR-2 process
    k = k+1;
  end
  act = mean(act,2);
  neur_act.soma = act*spike_opts.scalevec(1);
  neur_act.dend = act*spike_opts.scalevec(2);
  neur_act.bg = act*spike_opts.scalevec(3);
  [cleanIm,z_off] = scan_volume_frame(scan_vol,neur_act,scan_params,z_off);

  h5write(datapath,'/spikes',spikes,[1 k],[spike_opts.K numpts]);
  h5write(datapath,'/act',act,[1 i],[spike_opts.K 1]);
  h5write(datapath,'/motion',[0;0;z_off],[1 i],[3 1]);
  h5write(datapath,'/cleanimg',cleanIm,[1 1 i],[scan_params.vol_sz(1:2)/scan_params.sfrac 1]);
end

