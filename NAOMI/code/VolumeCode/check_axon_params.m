function axon_params = check_axon_params(axon_params)

% axon_params = check_axon_params(axon_params)
%  
% This function checks the elements of the struct axon_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - axon_params   - Struct containing parameters for axon generation
%       .flag        - Flag for generation of background dendrites
%       .distsc      - Parameter to determine how directed the random walk 
%                      to generate background processes is Higher values
%                      are more directed/less random (default = 0.5)
%       .fillweight  - Maximum length for a single process branch (default
%                      = 100 um)
%       .maxlength   - Maximum length for background processes (default =
%                      200 um) 
%       .minlength   - Minimum length for background processes (default =
%                      10 um) 
%       .maxdist     - Maximum distance to the end of a process (default =
%                      100 um) 
%       .maxel       - Max number of axons per voxel (default = 8)
%       .maxvoxel    - Max number of allowable axon elements per voxel
%                      (default = 6)
%       .numbranches - Number of allowable branches for a single process
%                      (default = 20) 
%       .varbranches - Standard deviation of the number of branches per
%                      process (default = 5) 
%       .maxfill     - Voxel maximum occupation: fraction of volume that is
%                      to be filled by background processes (default = 0.7)
%       .N_proc      - Number of background components (default = 10)
%       .l           - Gaussian process length scale for correllation
%                      background processes over the image (default = 25) 
%       .rho         - Gaussian process variance parameter for correllation
%                      background processes over the image (default = 0.1) 
%
% 2017 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(axon_params)                                                      % Make sure that axon_params is a struct
    clear axon_params
    axon_params = struct;
end

if (~isfield(axon_params,'flag'))||isempty(axon_params.flag)
    axon_params.flag = 1;
end
if (~isfield(axon_params,'distsc'))||isempty(axon_params.distsc)
    axon_params.distsc = 0.5;
end
if (~isfield(axon_params,'fillweight'))||isempty(axon_params.fillweight)
    axon_params.fillweight = 100;
end
if (~isfield(axon_params,'maxlength'))||isempty(axon_params.maxlength)
    axon_params.maxlength = 200;
end
if (~isfield(axon_params,'minlength'))||isempty(axon_params.minlength)
    axon_params.minlength = 10;
end
if (~isfield(axon_params,'maxdist'))||isempty(axon_params.maxdist)
    axon_params.maxdist = 100;
end
if (~isfield(axon_params,'maxel'))||isempty(axon_params.maxel)
    axon_params.maxel = 8;
end
if (~isfield(axon_params,'varfill'))||isempty(axon_params.varfill)             % variation in filling weight (std around 1)
    axon_params.varfill = 0.3;
end
if (~isfield(axon_params,'maxvoxel'))||isempty(axon_params.maxvoxel)           % Maximum number of elements per voxel
    axon_params.maxvoxel = 6;
end
if (~isfield(axon_params,'padsize'))||isempty(axon_params.padsize)             % Background padding size (for smoothness in background)
    axon_params.padsize = 20;
end
if (~isfield(axon_params,'numbranches'))||isempty(axon_params.numbranches)
    axon_params.numbranches = 20;
end
if (~isfield(axon_params,'varbranches'))||isempty(axon_params.varbranches)
    axon_params.varbranches = 5;
end
if (~isfield(axon_params,'maxfill'))||isempty(axon_params.maxfill)
    axon_params.maxfill = 0.5;
end
if (~isfield(axon_params,'N_proc'))||isempty(axon_params.N_proc)
    axon_params.N_proc = 10;
end
if (~isfield(axon_params,'l'))||isempty(axon_params.l)
    axon_params.l = 25;
end
if (~isfield(axon_params,'rho'))||isempty(axon_params.rho)
    axon_params.rho = 0.1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%