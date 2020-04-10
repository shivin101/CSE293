function bg_params = check_bg_params(bg_params)

% bg_params = check_bg_params(bg_params)
%  
% This function checks the elements of the struct bg_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - bg_params   - Struct containing parameters for background generation
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
%
% 2017 - Adam Charles and Alex Song
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(bg_params)                                                      % Make sure that bg_params is a struct
    clear bg_params
    bg_params = struct;
end

if (~isfield(bg_params,'flag'))||isempty(bg_params.flag)
    bg_params.flag = 1;
end
if (~isfield(bg_params,'distsc'))||isempty(bg_params.distsc)
    bg_params.distsc = 0.5;
end
if (~isfield(bg_params,'fillweight'))||isempty(bg_params.fillweight)
    bg_params.fillweight = 100;
end
if (~isfield(bg_params,'maxlength'))||isempty(bg_params.maxlength)
    bg_params.maxlength = 200;
end
if (~isfield(bg_params,'minlength'))||isempty(bg_params.minlength)
    bg_params.minlength = 10;
end
if (~isfield(bg_params,'maxdist'))||isempty(bg_params.maxdist)
    bg_params.maxdist = 100;
end
if (~isfield(bg_params,'maxel'))||isempty(bg_params.maxel)
    bg_params.maxel = 1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%