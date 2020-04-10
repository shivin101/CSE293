%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Plotting function for algorithmic comparisons %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('/media/adam/BigStorage/Calcium_Data/NAOMi_simulations/20180730_fullAlgoWorkspace.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if ~exist('ideals','var')
    ideals = struct();
end
cutoff = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get correlations for all the component traces

[cnmf, pcaica, suite2p, est] = correlateAllSegmentations(cnmf, pcaica, suite2p, est, neur_act2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make plot of all componants overlapping

[cnmf,suite2p,pcaica,est] = separateByCorrelation(cnmf,suite2p,pcaica,est,cutoff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make an image of all the ideal components separated by pair strength

ideals = separateIdealByCorrelation(ideals,est,cnmf,pcaica,suite2p,cutoff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find all the profiles that doubled up (explain the same ideal prof)

ideals = findDoubleProfiles(ideals, cnmf, suite2p, pcaica);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display some statistics

displayPairingStatistics(cnmf,pcaica,suite2p,est);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display all the strongly paired components (corr > 0.5)

comColors = [0.8,0.53,0;      % CNMF
             0.12,0.2,0.12;   % PCAICA
             0.2,0.5,0.6;     % Suite2p
             0.5,0.5,0.5;     % All components
             1,0.4,0.1;       % Ideal
             0.1,1,1];        % Ideal NOT found
         
% comColors = [1,0.57,0;    % CNMF
%              0.2,0.9,0.2; % PCAICA
%              0,0.3,1;     % Suite2p
%              0.5,0.5,0.5; % All components
%              1,1,0.1;     % 
%              0.1,1,1];    % Ideal NOT found

figure(1); 
subplot(3,6,[1,9]), imagesc((bsxfun(@times,sum(suite2p.allpairs.strongcomp,3), reshape(comColors(3,:),[1,1,3])) + ...
     bsxfun(@times,sum(cnmf.allpairs.strongcomp,3),reshape(comColors(1,:),[1,1,3])) + ...
     bsxfun(@times,sum(pcaica.allpairs.strongcomp,3),reshape(comColors(2,:),[1,1,3])))./max(1,...
     sum(suite2p.allpairs.strongcomp,3)>0 +  sum(cnmf.allpairs.strongcomp,3)>0 + ...
     sum(pcaica.allpairs.strongcomp,3)>0))
text(390,540, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)
axis image; axis off
subplot(3,6,[4,12]), imagesc(...
     bsxfun(@times,ideals.strongcomps, reshape(comColors(4,:),[1,1,3])) + ...
     bsxfun(@times,ideals.strongcompsI,reshape(comColors(5,:),[1,1,3])) + ...
     bsxfun(@times,ideals.strongcompsX,reshape(comColors(6,:),[1,1,3])))
text(390,540, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)
axis image; axis off
subplot(3,6,[13,14]), imagesc(bsxfun(@times,sum(cnmf.allpairs.strongcomp,3),reshape(comColors(1,:),[1,1,3])))
axis image; axis off
subplot(3,6,[15 16]), imagesc(bsxfun(@times,sum(suite2p.allpairs.strongcomp,3), reshape(comColors(3,:),[1,1,3])))
axis image; axis off
subplot(3,6,[17,18]), imagesc(bsxfun(@times,sum(pcaica.allpairs.strongcomp,3),reshape(comColors(2,:),[1,1,3])))
text(390, 570, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)
axis image; axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Show strong/weak/unpaired

Nr = 5;
Nc = 3;
figure(2)
subplot(Nr,Nc,1), imagesc(bsxfun(@times,ideals.strongcomps, reshape(comColors(4,:),[1,1,3])))
axis image; axis off
title('Strong pairings')
subplot(Nr,Nc,2), imagesc(bsxfun(@times,ideals.weakcomps, reshape(comColors(4,:),[1,1,3])))
axis image; axis off
title('Weak pairings')
subplot(Nr,Nc,3), imagesc(bsxfun(@times,ideals.npcomps, reshape(comColors(4,:),[1,1,3])))
axis image; axis off
title('Unpaired')
subplot(Nr,Nc,4), imagesc(bsxfun(@times,sum(cnmf.allpairs.strongcomp,3),reshape(comColors(1,:),[1,1,3]))) 
axis image; axis off
subplot(Nr,Nc,5), imagesc(bsxfun(@times,sum(cnmf.allpairs.weakcomp,3),reshape(comColors(1,:),[1,1,3])))
axis image; axis off
subplot(Nr,Nc,6), imagesc(bsxfun(@times,sum(cnmf.allpairs.upcomp,3),reshape(comColors(1,:),[1,1,3])))
axis image; axis off
subplot(Nr,Nc,7), imagesc(bsxfun(@times,sum(pcaica.allpairs.strongcomp,3),reshape(comColors(2,:),[1,1,3])))
axis image; axis off
subplot(Nr,Nc,8), imagesc(bsxfun(@times,sum(pcaica.allpairs.weakcomp,3),reshape(comColors(2,:),[1,1,3])))
axis image; axis off
subplot(Nr,Nc,9), imagesc(bsxfun(@times,sum(pcaica.allpairs.upcomp,3),reshape(comColors(2,:),[1,1,3])))
axis image; axis off
subplot(Nr,Nc,10), imagesc(bsxfun(@times,sum(suite2p.allpairs.strongcomp,3), reshape(comColors(3,:),[1,1,3])))
axis image; axis off
subplot(Nr,Nc,11), imagesc(bsxfun(@times,sum(suite2p.allpairs.weakcomp,3), reshape(comColors(3,:),[1,1,3])))
axis image; axis off
subplot(Nr,Nc,12), imagesc(bsxfun(@times,sum(suite2p.allpairs.upcomp,3), reshape(comColors(3,:),[1,1,3]))) 
axis image; axis off

subplot(Nr,Nc,13), imagesc(bsxfun(@times,est.allpairs.strongcomp,reshape(comColors(5,:),[1,1,3])))
axis image; axis off
subplot(Nr,Nc,14), imagesc(bsxfun(@times,est.allpairs.weakcomp,reshape(comColors(5,:),[1,1,3])))
axis image; axis off
subplot(Nr,Nc,15), imagesc(bsxfun(@times,est.allpairs.upcomp,reshape(comColors(5,:),[1,1,3])))
axis image; axis off
text(390, 570, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)

clear Nr Nc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Show boundary plots

figure(3); 
subplot(1,3,1), imagesc(0.5*repmat(ideals.strongcomps>0.2,[1,1,3]))
text(390,540, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)
hold on
for kk = 1:numel(cnmf.strongbound)
    if ~isempty(cnmf.strongbound{kk})
        plot(cnmf.strongbound{kk}(:,2),cnmf.strongbound{kk}(:,1),'r','LineWidth',2);
    end
end
hold off
axis image; axis off
subplot(1,3,2), imagesc(0.5*repmat(ideals.strongcomps>0.2,[1,1,3]))
text(390,540, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)
hold on
for kk = 1:numel(suite2p.strongbound)
    if ~isempty(suite2p.strongbound{kk})
        plot(suite2p.strongbound{kk}(:,2),suite2p.strongbound{kk}(:,1),'b','LineWidth',2);
    end
end
hold off
axis image; axis off
subplot(1,3,3), imagesc(0.5*repmat(ideals.strongcomps>0.2,[1,1,3]))
text(390,540, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)
hold on
for kk = 1:numel(pcaica.strongbound)
    if ~isempty(pcaica.strongbound{kk})
        plot(pcaica.strongbound{kk}(:,2),pcaica.strongbound{kk}(:,1),'g','LineWidth',2);
    end
end
hold off
axis image; axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot an image with doubling

subIm = subIndexImages(cnmf,suite2p, pcaica, est, ideals.doubling.cnmf(2));
% subIm = subIndexImages(cnmf,suite2p, pcaica, est, ideals.doubling.suite2p(4));
figure(3);
subplot(1,2,2), imagesc(0.5*subIm.cnmf+subIm.suite2p+0.5*subIm.pcaica)
text(390,540, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)
axis image; axis off
subplot(1,2,1), imagesc(2*repmat(subIm.ideals,[1,1,3]))
text(390,540, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)
hold on
for kk = 1:numel(subIm.cnmfbound)
    if ~isempty(subIm.cnmfbound{kk})
        plot(subIm.cnmfbound{kk}(:,2),subIm.cnmfbound{kk}(:,1),'r','LineWidth',2);
    end
end
for kk = 1:numel(subIm.suite2pbound)
    if ~isempty(subIm.suite2pbound{kk})
        plot(subIm.suite2pbound{kk}(:,2),subIm.suite2pbound{kk}(:,1),'b','LineWidth',2);
    end
end
for kk = 1:numel(subIm.pcaicabound)
    if ~isempty(subIm.pcaicabound{kk})
        plot(subIm.pcaicabound{kk}(:,2),subIm.pcaicabound{kk}(:,1),'g','LineWidth',2);
    end
end
hold off
axis image; axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot histograms of correlation values for all algorithms

figure(4);
edges = linspace(0.05, 1, 25);
h1 = histcounts(cnmf.corrvals,edges);
h2 = histcounts(suite2p.corrvals,edges);
h3 = histcounts(pcaica.corrvals,edges);
h4 = histcounts(est.corrvals,edges);
h  = bar(edges(1:end-1),[h4; h1; h2; h3]',1);
% comColors
set(h(1),'FaceColor',comColors(6,:)); 
set(h(2),'FaceColor',comColors(1,:)); 
set(h(3),'FaceColor',comColors(2,:));
set(h(4),'FaceColor',comColors(3,:));
% set(h(3),'FaceColor',get(h(1),'FaceColor'));
% set(h(4),'FaceColor',[0.3,0.7,0.1]);
% set(h(1),'FaceColor',[0.9,0.9,0.1]);
set(gca,'Xlim',[0.05,1])
box off
xlabel('Correlation of paried profiles','FontSize',18)
ylabel('Frequency','FontSize',18)
legend('Ideal','CNMF','Suite2p','PCA/ICA')

clear h1 h2 h3 h4 edges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROC curve

roc = genAlgorithmPairROC(cnmf,pcaica,suite2p,est,500,0.5);

figure(6)
plot(roc.FA.est,roc.TP.est,'k',roc.FA.cnmf,roc.TP.cnmf,'r',...
    roc.FA.suite2p ,roc.TP.suite2p,'b',roc.FA.pcaica ,roc.TP.pcaica,'g')
legend('Ideal', 'CNMF','Suite2p','PCA/ICA')
set(gca,'XLim',[0,1],'YLim',[0,1])
xlabel('False alarm rate')
ylabel('True positive rate')
axis square
box off
% subplot(1,2,2), plot(corr_thresh,roc.cnmf(:,1),'r',corr_thresh,roc.suite2p(:,1),'b',corr_thresh,roc.pcaica(:,1),'g')
% legend('CNMF','Suite2p','PCA/ICA')
% box off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot example time traces

% Fig 6: nsel: 1, 47
% nsel       = 1; t_ss       = 1000+(1:2*30*60);
nsel       = 50;
t_ss       = 1000+(1:2*30*60);
tt         = (0:size(est.estact,2)-1)/30;
all_strong = intersect(intersect(est.allpairs.strongpairs(:,1), cnmf.allpairs.strongpairs(:,1)),...
    intersect(suite2p.allpairs.strongpairs(:,1),pcaica.allpairs.strongpairs(:,1)));

figure(5);
subplot(5,1,1), plot(tt(t_ss),neur_act2(all_strong(nsel),t_ss)/max(neur_act2(all_strong(nsel),:)),'k')
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
box off
title('Ground truth')
subplot(5,1,2), plot(tt(t_ss),est.estact(est.allpairs.strongpairs(est.allpairs.strongpairs(:,1)==all_strong(nsel),2),t_ss)./ ...
                  max(est.estact(est.allpairs.strongpairs(est.allpairs.strongpairs(:,1)==all_strong(nsel),2),:)),'color',[0.9,0.9,0.1])
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
box off
title(sprintf('Ideal correlation: %1.4f',est.corrvals(est.allpairs.strongpairs(:,1)==all_strong(nsel))))
subplot(5,1,3), plot(tt(t_ss),cnmf.compTimecourse(cnmf.allpairs.strongpairs(cnmf.allpairs.strongpairs(:,1)==all_strong(nsel),2),t_ss)./...
    max(vec(cnmf.compTimecourse(cnmf.allpairs.strongpairs(cnmf.allpairs.strongpairs(:,1)==all_strong(nsel),2),:))),'r')
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title(sprintf('CNMF correlation: %1.4f',cnmf.corrvals(cnmf.allpairs.strongpairs(:,1)==all_strong(nsel))))
box off
subplot(5,1,4), plot(tt(t_ss),suite2p.compTimecourse(suite2p.allpairs.strongpairs(suite2p.allpairs.strongpairs(:,1)==all_strong(nsel),2),t_ss)./...
    max(vec(suite2p.compTimecourse(suite2p.allpairs.strongpairs(suite2p.allpairs.strongpairs(:,1)==all_strong(nsel),2),:))),'b')
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title(sprintf('Suite2p correlation: %1.4f',suite2p.corrvals(suite2p.allpairs.strongpairs(:,1)==all_strong(nsel))))
box off
subplot(5,1,5), plot(tt(t_ss),pcaica.compTimecourse(pcaica.allpairs.strongpairs(pcaica.allpairs.strongpairs(:,1)==all_strong(nsel),2),t_ss)./...
    max(vec(pcaica.compTimecourse(pcaica.allpairs.strongpairs(pcaica.allpairs.strongpairs(:,1)==all_strong(nsel),2),:))),'g')
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title(sprintf('PCA/ICA correlation: %1.4f',pcaica.corrvals(pcaica.allpairs.strongpairs(:,1)==all_strong(nsel))))
box off


clear all_strong tt nsel t_ss


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

m4x = matfile('/home/adam/GITrepos/tao_sim/Results/20180802_fullAlgoWorkspace_x4.mat');
TMPest     = m4x.est;
TMPcnmf    = m4x.cnmf;
TMPsuite2p = m4x.suite2p;
TMPpcaica  = m4x.pcaica;

% clear TMPest
% clear TMPcnmf
% clear TMPsuite2p
% clear TMPpcaica

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot example time traces across powers

% m4x = matfile('/home/adam/GITrepos/tao_sim/Results/20180802_fullAlgoWorkspace_x4.mat');
nsel       = 51;
t_ss       = 1000+(1:2*30*60);
tt         = (0:size(est.estact,2)-1)/30;
all_strong = intersect(intersect(est.allpairs.strongpairs(:,1), cnmf.allpairs.strongpairs(:,1)),...
    intersect(suite2p.allpairs.strongpairs(:,1),pcaica.allpairs.strongpairs(:,1)));

figure(5);
subplot(5,2,1), plot(tt(t_ss),neur_act2(all_strong(nsel),t_ss)/max(neur_act2(all_strong(nsel),:)),'k')
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title('Ground truth')
box off
subplot(5,2,3), plot(tt(t_ss),est.estact(est.allpairs.strongpairs(est.allpairs.strongpairs(:,1)==all_strong(nsel),2),t_ss)./ ...
                  max(est.estact(est.allpairs.strongpairs(est.allpairs.strongpairs(:,1)==all_strong(nsel),2),:)),'color',[0.9,0.9,0.1])
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title(sprintf('Ideal correlation: %1.4f',est.corrvals(est.allpairs.strongpairs(:,1)==all_strong(nsel))))
box off
subplot(5,2,5), plot(tt(t_ss),cnmf.compTimecourse(cnmf.allpairs.strongpairs(cnmf.allpairs.strongpairs(:,1)==all_strong(nsel),2),t_ss)./...
    max(vec(cnmf.compTimecourse(cnmf.allpairs.strongpairs(cnmf.allpairs.strongpairs(:,1)==all_strong(nsel),2),:))),'r')
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title(sprintf('CNMF correlation: %1.4f',cnmf.corrvals(cnmf.allpairs.strongpairs(:,1)==all_strong(nsel))))
box off
subplot(5,2,7), plot(tt(t_ss),suite2p.compTimecourse(suite2p.allpairs.strongpairs(suite2p.allpairs.strongpairs(:,1)==all_strong(nsel),2),t_ss)./...
    max(vec(suite2p.compTimecourse(suite2p.allpairs.strongpairs(suite2p.allpairs.strongpairs(:,1)==all_strong(nsel),2),:))),'b')
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title(sprintf('Suite2p correlation: %1.4f',suite2p.corrvals(suite2p.allpairs.strongpairs(:,1)==all_strong(nsel))))
box off
subplot(5,2,9), plot(tt(t_ss),pcaica.compTimecourse(pcaica.allpairs.strongpairs(pcaica.allpairs.strongpairs(:,1)==all_strong(nsel),2),t_ss)./...
    max(vec(pcaica.compTimecourse(pcaica.allpairs.strongpairs(pcaica.allpairs.strongpairs(:,1)==all_strong(nsel),2),:))),'g')
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title(sprintf('PCA/ICA correlation: %1.4f',pcaica.corrvals(pcaica.allpairs.strongpairs(:,1)==all_strong(nsel))))
box off


figure(5);
subplot(5,2,2), plot(tt(t_ss),neur_act2(all_strong(nsel),t_ss)/max(neur_act2(all_strong(nsel),:)),'k')
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title('Ground truth')
box off
% TMPest = m4x.est;
subplot(5,2,4), plot(tt(t_ss),TMPest.estact(TMPest.allpairs.strongpairs(TMPest.allpairs.strongpairs(:,1)==all_strong(nsel),2),t_ss)./ ...
                  max(TMPest.estact(TMPest.allpairs.strongpairs(TMPest.allpairs.strongpairs(:,1)==all_strong(nsel),2),:)),'color',[0.9,0.9,0.1])
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title(sprintf('Ideal correlation: %1.4f',TMPest.corrvals(TMPest.allpairs.strongpairs(:,1)==all_strong(nsel))))
box off
% clear TMPest
% TMPcnmf = m4x.cnmf;
subplot(5,2,6), plot(tt(t_ss),TMPcnmf.compTimecourse(TMPcnmf.allpairs.strongpairs(TMPcnmf.allpairs.strongpairs(:,1)==all_strong(nsel),2),t_ss)./...
    max(vec(TMPcnmf.compTimecourse(TMPcnmf.allpairs.strongpairs(TMPcnmf.allpairs.strongpairs(:,1)==all_strong(nsel),2),:))),'r')
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title(sprintf('CNMF correlation: %1.4f',TMPcnmf.corrvals(TMPcnmf.allpairs.strongpairs(:,1)==all_strong(nsel))))
box off
% clear TMPcnmf
% TMPsuite2p = m4x.suite2p;
subplot(5,2,8), plot(tt(t_ss),TMPsuite2p.compTimecourse(TMPsuite2p.allpairs.strongpairs(TMPsuite2p.allpairs.strongpairs(:,1)==all_strong(nsel),2),t_ss)./...
    max(vec(TMPsuite2p.compTimecourse(TMPsuite2p.allpairs.strongpairs(TMPsuite2p.allpairs.strongpairs(:,1)==all_strong(nsel),2),:))),'b')
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title(sprintf('Suite2p correlation: %1.4f',TMPsuite2p.corrvals(TMPsuite2p.allpairs.strongpairs(:,1)==all_strong(nsel))))
box off
% clear TMPsuite2p
% TMPpcaica = m4x.pcaica;
subplot(5,2,10), plot(tt(t_ss),TMPpcaica.compTimecourse(TMPpcaica.allpairs.strongpairs(TMPpcaica.allpairs.strongpairs(:,1)==all_strong(nsel),2),t_ss)./...
    max(vec(TMPpcaica.compTimecourse(TMPpcaica.allpairs.strongpairs(TMPpcaica.allpairs.strongpairs(:,1)==all_strong(nsel),2),:))),'g')
set(gca,'YTick',[],'XLim',[min(tt(t_ss)),max(tt(t_ss))])
title(sprintf('PCA/ICA correlation: %1.4f',TMPpcaica.corrvals(TMPpcaica.allpairs.strongpairs(:,1)==all_strong(nsel))))
box off
% clear TMPpcaica

clear all_strong tt nsel t_ss

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ideal Supp Fig

figure(3); 
subplot(1,3,1), imagesc(0.5*repmat(ideals.strongcomps>0.2,[1,1,3]))
text(390,540, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)
hold on
for kk = 1:numel(est.strongbound)
    if ~isempty(est.strongbound{kk})
        plot(est.strongbound{kk}(:,2),est.strongbound{kk}(:,1),'color',[0.9,0.9,0.1],'LineWidth',2);
    end
end
hold off
title('Paired profiles')
axis image; axis off

subplot(1,3,[2,3]), imagesc(est.estact(est.pairs(est.corrvals>cutoff,2),1000+(1:30*60)),[0,16])
axis image
colormap(fireprint)
set(gca,'XTick',[0,30*30,30*60],'XTickLabel',[0,30,60],'YDir','normal')
title('Paired time-traces')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CNMF Supp Fig


figure(3); 
subplot(1,3,1), imagesc(0.5*repmat(ideals.strongcomps>0.2,[1,1,3]))
text(390,540, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)
hold on
for kk = 1:numel(cnmf.strongbound)
    if ~isempty(cnmf.strongbound{kk})
        plot(cnmf.strongbound{kk}(:,2),cnmf.strongbound{kk}(:,1),'r','LineWidth',2);
    end
end
hold off
title('Paired profiles')
axis image; axis off

subplot(1,3,[2,3]), imagesc(cnmf.compTimecourse(cnmf.pairs(cnmf.corrvals>cutoff,2),1000+(1:30*60)),[0,10])
axis image
colormap(fireprint)
set(gca,'XTick',[0,30*30,30*60],'XTickLabel',[0,30,60],'YDir','normal')
title('Paired time-traces')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Suite2p Supp Fig

figure(3); 
subplot(1,4,1), imagesc(0.5*repmat(ideals.strongcomps>0.2,[1,1,3]))
text(390,540, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)
hold on
for kk = 1:numel(suite2p.strongbound)
    if ~isempty(suite2p.strongbound{kk})
        plot(suite2p.strongbound{kk}(:,2),suite2p.strongbound{kk}(:,1),'b','LineWidth',2);
    end
end
hold off
title('Paired profiles')
axis image; axis off

subplot(1,3,[2,3]), imagesc(suite2p.compTimecourse(suite2p.pairs(suite2p.corrvals>cutoff,2),1000+(1:30*60)),[0,700])
axis image
colormap(fireprint)
set(gca,'XTick',[0,30*30,30*60],'XTickLabel',[0,30,60],'YDir','normal')
title('Paired time-traces')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PCA/ICA Supp Fig

figure(3); 
subplot(1,4,1), imagesc(0.5*repmat(ideals.strongcomps>0.2,[1,1,3]))
text(390,540, '100 um','FontSize', 18)
line([390 490], 530*[1 1],'color','k','LineWidth',4)
hold on
for kk = 1:numel(pcaica.strongbound)
    if ~isempty(pcaica.strongbound{kk})
        plot(pcaica.strongbound{kk}(:,2),pcaica.strongbound{kk}(:,1),'g','LineWidth',2);
    end
end
hold off
title('Paired profiles')
axis image; axis off

subplot(1,3,[2,3]), imagesc(pcaica.compTimecourse(pcaica.pairs(pcaica.corrvals>cutoff,2),1000+(1:30*60)),[0,700])
axis image
colormap(fireprint)
set(gca,'XTick',[0,30*30,30*60],'XTickLabel',[0,30,60],'YDir','normal')
title('Paired time-traces')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% 
TMP.suite2pTT = bsxfun(@times, vec(sum(sum(suite2p.compSpatial>0))), bsxfun(@plus,  -median(suite2p.compTimecourse,2), suite2p.compTimecourse));

TMP.max       = max(max(TMP.suite2pTT));
TMP.allmaxs   = max(TMP.suite2pTT,[],2);

[~,TMP.IX] = sort(TMP.allmaxs,'descend');
TMP.IXlable = zeros(size(TMP.IX));
for ll = 1:numel(TMP.IX)
    if sum(TMP.IX(ll)==suite2p.allpairs.strongpairs) > 0
        TMP.IXlable(ll) = 1;
    else
        TMP.IXlable(ll) = 0;
    end
end

% clear TMP ll

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary figure (Zoomed segmentation)

figure(2);
% isub = [275,310,90,125];
% isub = [415,437,420,440];
% isub = [206, 256, 56, 106];
isub = [245, 293, 260, 305];
subplot(5,6,[1,15]), imagesc(bsxfun(@times,sum(ideals.strongcomps(isub(1):isub(2),isub(3):isub(4),:),3), ...
     reshape(comColors(4,:),[1,1,3])))
text(1, 11+(isub(4)-isub(3)), '10 um','FontSize', 18)
line([1 10], 8+(isub(4)-isub(3))*[1 1],'color','k','LineWidth',4)
axis image
axis off
subplot(5,6,[4,18]), imagesc(bsxfun(@times,sum(cnmf.allpairs.strongcomp(isub(1):isub(2),isub(3):isub(4),:),3),reshape(comColors(1,:),[1,1,3]))...
             + bsxfun(@times,sum(suite2p.allpairs.strongcomp(isub(1):isub(2),isub(3):isub(4),:),3),reshape(comColors(3,:),[1,1,3]))...
             + bsxfun(@times,sum(pcaica.allpairs.strongcomp(isub(1):isub(2),isub(3):isub(4),:),3),reshape(comColors(2,:),[1,1,3])))
text(1, 11+(isub(4)-isub(3)), '10 um','FontSize', 18)
line([1 10], 8+(isub(4)-isub(3))*[1 1],'color','k','LineWidth',4)
axis image
axis off
subplot(5,6,[19,26]), imagesc(bsxfun(@times,sum(cnmf.allpairs.strongcomp(isub(1):isub(2),isub(3):isub(4),:),3),reshape(comColors(1,:),[1,1,3])))
axis image
axis off
subplot(5,6,[21,28]), imagesc(bsxfun(@times,sum(suite2p.allpairs.strongcomp(isub(1):isub(2),isub(3):isub(4),:),3),reshape(comColors(3,:),[1,1,3])))
axis image
axis off
subplot(5,6,[23,30]), imagesc(bsxfun(@times,sum(pcaica.allpairs.strongcomp(isub(1):isub(2),isub(3):isub(4),:),3),reshape(comColors(2,:),[1,1,3])))
text(1, 11+(isub(4)-isub(3)), '10 um','FontSize', 18)
line([1 10], 8+(isub(4)-isub(3))*[1 1],'color','k','LineWidth',4)
axis image
axis off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
