function tpm_params = check_tpm_params(tpm_params)

% tpm_params = check_tpm_params(tpm_params)
%  
% This function checks the elements of the struct vol_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - tpm_params - Struct containing parameters for estimating photon flux
%     .nidx   = 1.33       - refractive index [], water
%     .nac    = 0.8        - objective NA []
%     .phi    = 0.8*SA*0.4 - 80 percent transmission and 0.8NA solid angle
%                            with average 40 perfect PMT QE []
%     .eta    = 0.6        - eGFP quantum yield, estimate for 2P QY []
%     .conc   = 10         - fluorophore concentration, average of  
%                            literature (Huber et al 2012, Zariwala et al
%                            2012), [uM] 
%     .delta  = 35         - two-photon abs. cross section, estimate from
%                            Harris lab (Janelia), saturated GCaMP [GM]
%     .delta  = 2          - estimate of 2pCS at resting calcium levels
%                            [GM] 
%     .gp     = 0.588      - pulse-shape temporal coherence [], sech 
%     .f      = 80         - Ti:S laser rep rate [MHz]
%     .tau    = 150        - Ti:S laser pulse-width [fs]
%     .pavg   = 40         - laser power [mW]
%     .lambda = 0.92       - wavelength used for GCaMP excitation
% 
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(tpm_params)                                                     % Make sure that tpm_params is a struct
    clear tpm_params
    tpm_params = struct;
end

if (~isfield(tpm_params,'nidx'))||isempty(tpm_params.nidx)                 % water
    tpm_params.nidx = 1.33;                                            
end

if (~isfield(tpm_params,'nac'))||isempty(tpm_params.nac)                   % objective NA
    tpm_params.nac = 0.8;                                            
end

if (~isfield(tpm_params,'phi'))||isempty(tpm_params.phi)                   % 80 percent transmission and 0.8NA light all on detector (no strong scattering) with average 40 perfect PMT QE
    tpm_params.phi = 0.8*((1-sqrt(1-(tpm_params.nac/tpm_params.nidx)^2))/2)*0.4;                                           
end

if (~isfield(tpm_params,'eta'))||isempty(tpm_params.eta)                   % eGFP quantum yield
    tpm_params.eta = 0.6;                                            
end

if (~isfield(tpm_params,'conc'))||isempty(tpm_params.conc)                 %average of literature measurements (Huber et al 2012, Zariwala et al 2012), uM
    tpm_params.conc = 10;                                            
end

if (~isfield(tpm_params,'delta'))||isempty(tpm_params.delta)               %estimate from Harris lab (Janelia), saturated GCaMP is 35, estimate of F0 at resting calcium levels
    tpm_params.delta = 2;                                            
end

if (~isfield(tpm_params,'gp'))||isempty(tpm_params.gp)                     % sech pulse shape (temporal)
    tpm_params.gp = 0.588;                                            
end

if (~isfield(tpm_params,'f'))||isempty(tpm_params.f)                       % Ti:S laser
    tpm_params.f = 80;                                            
end

if (~isfield(tpm_params,'tau'))||isempty(tpm_params.tau)                   % Ti:S laser
    tpm_params.tau = 150;                                            
end

if (~isfield(tpm_params,'pavg'))||isempty(tpm_params.pavg)                 % 50mW power for L2/3 imaging
    tpm_params.pavg = 40;                                            
end

if (~isfield(tpm_params,'lambda'))||isempty(tpm_params.lambda)             % 0.92um for excitation wavelength
    tpm_params.lambda = 0.92;                                            
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%