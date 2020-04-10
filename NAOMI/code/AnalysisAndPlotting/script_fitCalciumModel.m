%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit the NAOMi calcium model
% 
% Script that fits the NAOMi calcium model to simultaneously recorded
% spikes and fluorescence. The data used here was originally published in 
% (Chen et. al. 2013 Nature; Akerboom, Chen 2012 J. Neurosci). This script
% was in part adapted from a script written by Tsai-Wen Chen (2015/01/27)
% that called the data from the above publications. 
%
% From the original script:
% ------------------------------------------
% Ephys data were recorded at 10KHz
% Imaging data were recorded at 60Hz
%
% Each .mat data file contains a variable named 'obj'
% key recording traces and time base can be accessed by the following:
%
% traces = obj.timeSeriesArrayHash.value{id}.valueMatrix
% time   = obj.timeSeriesArrayHash.value{id}.time
%
% id=
% 
% 1: fmean_roi
% 2: fmean_neuropil
% 3: raw_ephys
% 4: filtered_ephys
% 5: detected_spikes
% ------------------------------------------
%
% 
% Adam Charles - 2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data

clear
load('/home/adam/GITrepos/tao_sim/Data/chenData/processed_data/data_20120521_cell2_003.mat') % Change this path to where the data is

fmean_roi       = obj.timeSeriesArrayHash.value{1}.valueMatrix;
fmean_neuropil  = obj.timeSeriesArrayHash.value{2}.valueMatrix;
fmean_comp      = fmean_roi-0.7*fmean_neuropil;                            % their estimate of fluorescence from cell
t_frame         = obj.timeSeriesArrayHash.value{1}.time;
filt            = obj.timeSeriesArrayHash.value{4}.valueMatrix;
t_ephys         = obj.timeSeriesArrayHash.value{4}.time;
detected_spikes = obj.timeSeriesArrayHash.value{5}.valueMatrix;
spike_time      = t_ephys(detected_spikes);
fmean_norm      = fmean_comp./quantile(fmean_comp,0.2);                    % normalized fluoresence to baseline of 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pick a subsample of the data

fmean_normCut = fmean_norm(1:3600);
spike_timeCut = spike_time(spike_time<=60);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run fitting and simulate the time-traces

stimes          = ceil(spike_time*100);                                    % 100Hz estimate of spike times
spikes2         = zeros(1,round(max(t_ephys)));
spikes2(stimes) = 7.6e-6;
spike_opts      = check_spike_opts([]);

cal_params = fit_NAOMi_calcium(spike_timeCut, ...
                    t_frame(1:length(fmean_normCut)), fmean_normCut, ...
                        max(t_frame(1:length(fmean_normCut))), 'fmincon');
cal_params.sat_type = 'Ca_DE';                                             % Set the calcium dynamics to 'single' mode (simpler dynamics)
cal_params.dt       = 1/100;                                               % Set the calcium dynamics framerate

[~,~,simFluor]  = calcium_dynamics(spikes2, cal_params, 'GCaMP6');
simFluor         = simFluor/min(simFluor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot traces

figure(1);
cla
plot(t_frame,fmean_norm/median(fmean_norm), 'b', (1:length(simFluor))/100,simFluor, '.-r', (1:length(simFluor7))/100,simFluor7, '--k', 'LineWidth', 2)
% plot(t_frame(1:length(fmean_normCut)),fmean_normCut/median(fmean_normCut), 'b', (1:length(simFluor))/100,simFluor, '.-r', (1:length(simFluor7))/100,simFluor7, '--k', 'LineWidth', 2)
hold on;
stem(spike_time,7*ones(length(spike_time),1),'.k');
set(gca, 'FontSize', 18, 'TickDir', 'out', 'XLim', [0, length(fmean_norm)/60], 'YLim', [0.8, 6.9], 'YTick', [])
box off
xlabel('Time (s)','FontSize', 18)
ylabel('\Delta F/F (AU)','FontSize', 18)

set(gcf,'color',[1,1,1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%