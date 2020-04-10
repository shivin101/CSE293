function neur_params = check_neur_params(neur_params)

% neur_params = check_neur_params(neur_params)
%  
% This function checks the elements of the struct neur_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - neur_params - Struct containing parameters for neuron generation
%       .n_samps     - Number of sphere samples to use in creating the mesh
%                      for generating soma and nucleus shapes (default =
%                      200) 
%       .l_scale     - length-scale for the isotropic GP of the soma
%                      shapes. This controls the shape `bumpiness' (default
%                      = 90) 
%       .p_scale     - Overall variance of the isotropic GP of the soma 
%                      shape. (default = 90) 
%       .avg_rad     - Average radius of each neuron in um (default =
%                      6.6 um) 
%       .nuc_fluorsc - Potential fluorescence in the nucleus (default =
%                      0)
%       .min_thic    - Minimum cytoplasmic thickness (default = 0.8)
%       .eccen       - Maximum eccentricity of neuron (default = 0.35)
%       .exts        - Parameters dictating the max/min of the soma radii
%                      (Default = [0.75,1.7])
%       .nexts       - Parameters dictating the extent to shrink and smooth
%                      the nucleus (Default = [0.5, 1])
%       .neur_type   - Option for neuron type (Default 'pyr')
%
% 2017 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(neur_params)                                                    % Make sure that neur_params is a struct
    clear neur_params
    neur_params = struct;
end

if (~isfield(neur_params,'n_samps'))||isempty(neur_params.n_samps)         % Default numbers of samples on the sphere is 1000             
    neur_params.n_samps = 200;
end
if (~isfield(neur_params,'l_scale'))||isempty(neur_params.l_scale)         % Default length-scale is 105     
    neur_params.l_scale = 90;
end
if (~isfield(neur_params,'p_scale'))||isempty(neur_params.p_scale)         % Default correlation scaling is 95
    neur_params.p_scale = 1000;
%     neur_params.p_scale = 90;
end
if (~isfield(neur_params,'avg_rad'))||isempty(neur_params.avg_rad)         % Default average radius is 6.6      
    neur_params.avg_rad = 5.9;
%     neur_params.avg_rad = 6.8;
end
if (~isfield(neur_params,'nuc_rad'))||isempty(neur_params.nuc_rad)         % Default nuclear radius    
    neur_params.nuc_rad = [5.65 2.5];
end
if (~isfield(neur_params,'max_ang'))||isempty(neur_params.max_ang)         % Default maximum angle tilt    
    neur_params.max_ang = 20;
end

if (~isfield(neur_params,'plot_opt'))||isempty(neur_params.plot_opt)       % Default plotting setting is to plot (Default FALSE)
    neur_params.plot_opt = false;
end
if (~isfield(neur_params,'dendrite_tau'))||isempty(neur_params.dendrite_tau)% Default axon option (Default FALSE)
    neur_params.dendrite_tau = 50;
end
if (~isfield(neur_params,'nuc_fluorsc'))||isempty(neur_params.nuc_fluorsc) % Nuclear fluoresence level
    neur_params.nuc_fluorsc = 0;
end
if (~isfield(neur_params,'min_thic'))||isempty(neur_params.min_thic)       % Minimum cytoplasmic thickness 
    neur_params.min_thic = [0.4 0.4];
%     neur_params.min_thic = 0.8;
end
if (~isfield(neur_params,'eccen'))||isempty(neur_params.eccen)             % Default maximum eccentricity of neuron is 0.25
    neur_params.eccen = [0.35 0.35 0.5];
%     neur_params.eccen = 0.35;
end
if (~isfield(neur_params,'exts'))||isempty(neur_params.exts)               % Parameters dictating the max/min of the soma radii
    neur_params.exts = [0.75,1.7];
end
if (~isfield(neur_params,'nexts'))||isempty(neur_params.nexts)             % Parameters dictating the extent to shrink and smooth a nucleus
    neur_params.nexts = [0.5,1];
end
if (~isfield(neur_params,'neur_type'))||isempty(neur_params.neur_type)     % Option for neuron type (Default 'pyr')
    neur_params.neur_type = 'pyr';
end
if (~isfield(neur_params,'fluor_dist'))||isempty(neur_params.fluor_dist)   % Somatic neural fluoresence distribution (mean length [um], mean weight [0-1])
%     neur_params.fluor_dist = [0 0];
    neur_params.fluor_dist = [1 0.2];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%