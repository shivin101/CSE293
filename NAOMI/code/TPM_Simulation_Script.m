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
% contacting the creators, Adam Charles (adamsc@princeton.edu) and Alex
% Song. 
% 
% 2019 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all paths
% The following script adds the paths to the various NAOMi scripts and
% functions. Furthermore it will run the mex_compiling script if any of the
% mex files needed were not already compiled. You can comment this out if
% you prefer to add the files and mex the files separately. 

% installNAOMi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Parameters
% This is an example of setting the parameters for simulated imaging. This
% simulation will generate a 100x100x100 micron volume 100 microns deep. It
% will then scan with an 0.6NA Gaussian beam using an 0.8NA objective with
% 40mW power for 5000 frames at 30Hz. It will then save the video to
% 'TMP_volume.mat'.

save_filename        = 'TPM_volume.tif';                                   % Pick a file to save the simulated movie as (ends with either .fits, .tif, or .mat)
vol_params.vol_sz    = [500,500,100];                                      % Volume size to sample (in microns)
vol_params.vol_depth = 100;                                                % Set the depth of imaging (depth at the middle of the simulated volume)
psf_params.objNA     = 0.8;                                                % Numerical aperture of PSF
psf_params.NA        = 0.6;                                                % Numerical aperture of PSF
tpm_params.pavg      = 40;                                                 % power in units of mW
spike_opts.nt        = 500;                                              % Set number of time step
spike_opts.dt        = 1/30;                                               % Sampling frame-rate period (1/Hz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check all the parameter structs to make sure all fields are set
% This block is not strickly needed as these checks are all run inside the
% necessary functions, however it may be informative to the user to see
% what the default parameters for each struct are set to 

vol_params   = check_vol_params(vol_params);                               % Check volume parameters
vasc_params  = check_vasc_params([]);                                      % Make default set of vasculature parameters
neur_params  = check_neur_params([]);                                      % Make default set of neuron parameters
dend_params  = check_dend_params([]);                                      % Make default set of dendrite parameters
axon_params  = check_axon_params([]);                                      % Make default set of axon parameters
bg_params    = check_bg_params([]);                                        % Make default set of background parameters
spike_opts   = check_spike_opts(spike_opts);                               % Check spike/fluorescence simulation parameters
noise_params = check_noise_params([]);                                     % Make default noise parameter struct for missing elements
psf_params   = check_psf_params(psf_params);                               % Check point spread function parameters
scan_params  = check_scan_params([]);                                      % Check the scanning parameter struct
tpm_params   = check_tpm_params(tpm_params);                               % Check the auxiliary two-photon imaging parameter struct 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate neural volume
% The following function generates the anatomical volume to be imaged. The
% main output is vol_out, which contains a series of cell arrays that
% contain the locations and concentration of fluorescence for all active
% components in the volume. This struct is passed into the other functional
% blocks of the simulator and is the first piece of the simulation that
% should be run. 

tic
[vol_out,vol_params,neur_params,vasc_params,dend_params,bg_params, ...
    axon_params] = simulate_neural_volume(vol_params, neur_params, ...
            vasc_params, dend_params, bg_params, axon_params, psf_params); % Draw a random volume - this takes the longest amound of time
fprintf('Simulated neural volume in %f seconds.\n', toc); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate optical mask and PSF
% The following function simulates the point-spread function that will be
% used to scan the volume. The volume struct is used to generate a mask
% that accounts for how light popagates through the tissue and vasculature.

tic
PSF_struct = simulate_optical_propagation(vol_params,psf_params,vol_out);  % Create the point-spread function and mask for scanning
fprintf('Simulated optical propagation in %f seconds.\n', toc); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate temporal activity
% The following script creates fluorescence time-traces for each active
% component in the volume. These time-traces are generated using the same
% spiking, but slightly different parameters, for each of the somatic,
% dendritic and axon components. The resulting struct, neur_act, contains
% these time traces. The distances between neural locations is used to
% correlate the activity in a Hawkes model (recurrent network model similar
% to GLM). If desired, neur_act can be generated by taking time traces from
% other simulators or real data. 

tic
spike_opts.K      = size(vol_out.gp_vals,1);                               % Read off the number of neurons
spike_opts.rate   = 0.25;                                                  % The average rate of firing for all the components can be modulated using this parameter
[neur_act,spikes] = generateTimeTraces(spike_opts,[],vol_out.locs);        % Generate time traces using AR-2 process
fprintf('Simulated temporal activity in %f seconds.\n', toc); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform scanning
% The following function takes the outputs from all the previous functions
% and simulates the scanning procedure. The output contains both the noisy
% (Poisson-Gauss model or Dynode model and per-pixel bleedthrough)
% fluorescence movie and, if desired, a clean version of the movie that has
% no noise (but still has motion). Motion can be toggled on and off, for
% example, with the demonstrated parameter. 

tic
scan_params.motion = false;                                                % Toggle this parameter to turn motion simulation on/off
[Fsim,Fsim_clean]  = scan_volume(vol_out, PSF_struct, neur_act, ...
                       scan_params, noise_params, spike_opts, tpm_params); % Perform the scanning simulation
fprintf('Simulated scanning in %f seconds.\n', toc); 
                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Create "ground truth" spatial profiles and SNR-adjusted "ideal profiles"
% To calculate the "ground truth" spatial profiles, we zero out all neurons
% except for that individual neuron (the equivalent of only labeling one 
% neuron in the volume). We then scan the volume using the same procedure, 
% and repeat for all neurons.

[comps,baseim,ideal] = calculateIdealComps(vol_out,PSF_struct,neur_act,...
                                   scan_params,noise_params,spike_opts,...
                                                tpm_params, spike_opts.K); % Scan each neuron to make the ideal components
idealTraces          = times_from_profs(Fsim, ideal, -1);                  % Get the ideal time-traces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Save data to file
% There are a number of ways to save the resulting simulated data. The
% results can all be saved to a .mat file using the standard save function.
% Alternatively, the video can be exported to .tif or .fits files that can
% be read in and analyzed using other programming environments (i.e.,
% python or julia)

save C:\Users\shivi\OneDrive\Documents\Courses\CSE293\New_data\naomi_sim\20190729_LargeVolume_nomotion.mat -v7.3          % Saves all variables to a .mat file
write_TPM_movie(Fsim, save_filename);                                      % Saves Fsim to save_name which can be a .mat, .tif or .fits file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make a movie of the results
% The output can also be saved as an avi movie for display purposes. As an
% example, the following line of code saves the video to such a file.

make_avi(Fsim, './TMPVID.avi', 0.2);                                       % Make an avi of the noisy video
make_avi(Fsim_clean, './TMPVIDclean.avi', 0.2);                            % Make an avi of the clean video

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

