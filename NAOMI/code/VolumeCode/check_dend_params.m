function dend_params = check_dend_params(dend_params)

% dend_params = check_dend_params(dend_params)
%  
% This function checks the elements of the struct dend_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - dend_params - Struct containing parameters for dendrite simulation
%       .dtParams        - dendritic tree number,radius in um of branches
%                          (uniform across circle),radius in um of branches
%                          (uniform in z) (default = [40 150 50 1 10])
%       .atParams        - Apical dendrite number,radius in um of branches
%                          (uniform across circle),radius in um of branches
%                          (uniform in z),offset from center in um (default
%                          = [1 5 5 5 12]) 
%       .atParams2       - Through-volume apical dendrite number,radius in
%                          um of branches (uniform across circle),radius in
%                          um of branches (uniform in z),offset from center
%                          in um (default = = [1 5 5 5 18])
%       .dweight         - Weight for path planning randomness in the
%                          dendrites (default = 10) 
%       .bweight         - Weight for obstruction (default = 50)
%       .thicknessScale  - Scaling for dendrite thickness in um^2 (default
%                          = 0.5)
%       .dims            - dims set at 10 um/space (default = [30 30 30])
%       .dimsSS          - dims subsampling factor (10 samples per dim
%                          grid) (default = [10 10 10]) 
%       .rallexp         - Rall exponent that controls the cahnge in size
%                          of dendrites at branching locations (default =
%                          1.5) 
%       .weightScale     - scaling for dendrite fluorescence (1/distance
%                          scaling [um], distance scaling weight [0-1],
%                          variation in fluoresence (uniform) [0-1])
%                          (default = [150, 0.5, 0.5])
%
% 2017 - Adam Charles and Alex Song 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks 

if isempty(dend_params)                                                    % Make sure that dend_params is a struct
    clear dend_params
    dend_params = struct;
end

if (~isfield(dend_params,'dtParams'))||isempty(dend_params.dtParams)
    dend_params.dtParams       = [40 150 50 1 10];                         % dendritic tree number,radius in um of branches (uniform across circle),radius in um of branches (uniform in z)
end
if (~isfield(dend_params,'atParams'))||isempty(dend_params.atParams)
    dend_params.atParams       = [6 5 5 5 1];                              % apical dendrite (L2/3) number,radius in um of branches (uniform across circle),radius in um of branches (uniform in z),offset from center in um
%     dend_params.atParams       = [1 5 5 5 6];                              % apical dendrite (L2/3) number,radius in um of branches (uniform across circle),radius in um of branches (uniform in z),offset from center in um
end
if (~isfield(dend_params,'atParams2'))||isempty(dend_params.atParams)
    dend_params.atParams2      = [1 5 5 5 4];                              % apical dendrite (L5) number,radius in um of branches (uniform across circle),radius in um of branches (uniform in z),offset from center in um
%     dend_params.atParams2      = [1 5 5 5 9];                              % apical dendrite (L5) number,radius in um of branches (uniform across circle),radius in um of branches (uniform in z),offset from center in um
end
if (~isfield(dend_params,'dweight'))||isempty(dend_params.dweight)
    dend_params.dweight        = 10;                                       % weight for path planning randomness
end
if (~isfield(dend_params,'bweight'))||isempty(dend_params.bweight)
    dend_params.bweight        = 5;                                        % weight for obstruction
end
if (~isfield(dend_params,'thicknessScale'))||isempty(dend_params.thicknessScale)
    dend_params.thicknessScale = 0.5;                                      % scaling for dendrite thickness (int). Should be 1 for 1um sampling,(4 for 0.5um sampling)
end
if (~isfield(dend_params,'weightScale'))||isempty(dend_params.weightScale)
%     dend_params.weightScale = [150 0 0];                               % scaling for dendrite fluorescence (1/distance scaling [um], distance scaling weight [0-1], variation in fluoresence (uniform) [0-1])
    dend_params.weightScale = [150 1 0.8];                               % scaling for dendrite fluorescence (1/distance scaling [um], distance scaling weight [0-1], variation in fluoresence (uniform) [0-1])
end
if (~isfield(dend_params,'dims'))||isempty(dend_params.dims)
    dend_params.dims           = [60 60 60];                               % dims set at 10 um per space
end
if (~isfield(dend_params,'dimsSS'))||isempty(dend_params.dimsSS)
    dend_params.dimsSS         = [5 5 5];                                  % dims subsampling factor (10 samples per dim grid)
end
if (~isfield(dend_params,'rallexp'))||isempty(dend_params.rallexp)
    dend_params.rallexp         = 1.5;   
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%