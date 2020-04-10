function [cnmf, pcaica, suite2p, est] = correlateAllSegmentations(cnmf, pcaica, suite2p, est, neur_act2, varargin)

% [cnmf, pcaica, suite2p, est] = correlateAllSegmentations(cnmf, pcaica, suite2p, est, neur_act2, varargin)
%
% Function to pair the components from different segmentation algorithms
%
% 2018 - Adam Charles & Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 5
    calc_only = varargin{1};
else
    calc_only = [true, true, true, true];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CNMF

if calc_only(1)
    [cnmf.corrvals, cnmf.pairs, cnmf.corrM] = corrTimeAndSpace(...
                                     neur_act2',cnmf.compFluoresence',...
                                     est.compsIdealAll, cnmf.compSpatial); % Calculate pairs based on overlap and best temporal correlation
    cnmf.simTraces    = neur_act2(cnmf.pairs(:,1),:);                      % Extract the paired traces (simulated)
    cnmf.pairedTraces = cnmf.compFluoresence(cnmf.pairs(:,2),:);           % Extract the paired traces (estimated)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PCA/ICA

if calc_only(2)
    [pcaica.corrvals, pcaica.pairs, pcaica.corrM] = corrTimeAndSpace(...
                                   neur_act2',pcaica.compTimecourse',...
                                 est.compsIdealAll, pcaica.compSpatialSc); % Calculate pairs based on overlap and best temporal correlation
    pcaica.simTraces    = neur_act2(pcaica.pairs(:,1),:);                  % Extract the paired traces (simulated)
    pcaica.pairedTraces = pcaica.compTimecourse(pcaica.pairs(:,2),:);      % Extract the paired traces (estimated)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Suite2p

if calc_only(3)
    [suite2p.corrvals, suite2p.pairs, suite2p.corrM] = corrTimeAndSpace(...
                                 neur_act2', suite2p.compTimecourse', ... 
                                  est.compsIdealAll, suite2p.compSpatial); % Calculate pairs based on overlap and best temporal correlation
    suite2p.simTraces    = neur_act2(suite2p.pairs(:,1),:);                % Extract the paired traces (simulated)
    suite2p.pairedTraces = suite2p.compTimecourse(suite2p.pairs(:,2),:);   % Extract the paired traces (estimated)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ideal

if calc_only(4)
    [est.corrvals, est.pairs, est.corrM] = corrTimeAndSpace(...
                                neur_act2', est.estact(1:end-1,:)', ....
                                       est.compsIdealAll, est.compsIdeal); % Calculate pairs based on overlap and best temporal correlation
    est.simTraces    = neur_act2(est.pairs(:,1),:);                        % Extract the paired traces (simulated)
    est.pairedTraces = est.estact(est.pairs(:,2),:);                       % Extract the paired traces (estimated)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OLDER CODE
% 
% [cnmf.corrvals, cnmf.pairs, cnmf.corrM] = corrTimeAndSpace(...
%                          est.estact(1:end-1,:)',cnmf.compFluoresence',...
%                                      est.compsIdealAll, cnmf.compSpatial); % 
% cnmf.simTraces           = neur_act2(cnmf.pairs(:,1),:);
% cnmf.pairedTraces        = cnmf.compFluoresence(cnmf.pairs(:,2),:);
% 
% [pcaica.corrvals, pcaica.pairs, pcaica.corrM] = corrTimeAndSpace(...
%                         est.estact(1:end-1,:)',pcaica.compTimecourse',...
%                                  est.compsIdealAll, pcaica.compSpatialSc); % 
% pcaica.simTraces    = neur_act2(pcaica.pairs(:,1),:);
% pcaica.pairedTraces = pcaica.compTimecourse(pcaica.pairs(:,2),:);
% 
% [suite2p.corrvals, suite2p.pairs, suite2p.corrM] = corrTimeAndSpace(...
%                      est.estact(1:end-1,:)', suite2p.compTimecourse', ... 
%                                   est.compsIdealAll, suite2p.compSpatial); % 
% suite2p.simTraces    = neur_act2(suite2p.pairs(:,1),:);
% suite2p.pairedTraces = suite2p.compTimecourse(suite2p.pairs(:,2),:);





% [cnmf.corrvals, cnmf.pairs, cnmf.corrM] = corrTimeAndSpace(est.estact(1:end-1,:)',cnmf.compFluoresence',est.compsIdealAll, cnmf.compSpatial);
% cnmf.simTraces           = neur_act2(cnmf.pairs(:,1),:);
% cnmf.pairedTraces        = cnmf.compFluoresence(cnmf.pairs(:,2),:);
% 
% [pcaica.corrvals, pcaica.pairs, pcaica.corrM] = corrTimeAndSpace(est.estact(1:end-1,:)',pcaica.compTimecourse',est.compsIdealAll, pcaica.compSpatialSc);
% pcaica.simTraces    = neur_act2(pcaica.pairs(:,1),:);
% pcaica.pairedTraces = pcaica.compTimecourse(pcaica.pairs(:,2),:);
% 
% [suite2p.corrvals, suite2p.pairs, suite2p.corrM] = corrTimeAndSpace(est.estact(1:end-1,:)',suite2p.compTimecourse',est.compsIdealAll, suite2p.compSpatial);
% suite2p.simTraces    = neur_act2(suite2p.pairs(:,1),:);
% suite2p.pairedTraces = suite2p.compTimecourse(suite2p.pairs(:,2),:);