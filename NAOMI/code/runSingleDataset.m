%rootPath = 'Bezos-center/OB4_LargeFOVmicroscope/SimulationData';
% rootPath = '/mnt/Bezos-center/OB4_LargeFOVmicroscope/SimulationData';
rootPath = './';

% opts  = {'bessel'};                                                     % Set the different options
% conds = {'standard'};                                                   % Set the different conditions to run
opts  = {'bessel'};                                                       % Set the different options     
% conds = {'gcamp7','sparse','nuclear'};                                    % Set the different conditions to run      
% opts  = {'standard','bessel','stefo','vtwins'};                         % Set the different options     
conds = {'standard','gcamp7','sparse','nuclear'};                       % Set the different conditions to run      
filename = cell(length(opts)*length(conds),1);                            % Initialize the set of filenames 
kk = 0;                                                                   % Initialize filename counter
for ii = 1:length(opts)
    for jj = 1:length(conds)
        kk = kk+1;
        filename{kk} = ['params_' opts{ii} '_' conds{jj} '.mat'];
        TPM_Simulation_Parameters(fullfile(rootPath,filename{kk}), opts{ii}, conds{jj})
    end
end

%%
if ~exist('spikes','var')
    load(fullfile(rootPath,'20180504_spikes5000frames.mat'));
end
if ~exist('vol_out','var')
    load(fullfile(rootPath,'20180504_largeVolume.mat'));
end

kk = 0;
for ii = 1:length(opts)
    for jj = 1:length(conds)
        kk = kk+1;
        load(fullfile(rootPath,filename{kk}),'psf_params','spike_opts','scan_params','noise_params')
        load(fullfile(rootPath,filename{kk}),'PSF_struct')
        if ~exist('PSF_struct','var')
            PSF_struct = simulate_optical_propagation(vol_params,psf_params,vol_out);  % Create the point-spread function and mask for scanning
        end    
        neur_act                 = generateTimeTraces(spike_opts,spikes.somas);                   % Generate time traces using AR-2 process
        neur_act.bg              = neur_act.dend;
        [Fsim,Fsim_clean,motion] = scan_volume(vol_out, PSF_struct, ...
                                     neur_act, scan_params, noise_params); % Perform the scanning simulation
        save(fullfile(rootPath,['movie' filename{kk}(7:end)]),'-v7.3','PSF_struct','neur_act','Fsim','Fsim_clean','motion');
    end
end


