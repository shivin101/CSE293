clc;
clear all;
close all;

tic

ids=[2,8,4];

root = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\M1-2\';
ROI ={1,size(ids,2)};

%Hard coding some naming conventions in this stage of the process
for i = 1:size(ids,2)
    file(i) = 'KA0013_1801'+string(num2str(ids(i)+5,'%02d'))+'_001_001_summed_50.roi';
    
    %Create an array of ROI and images from preprocessed arrays
    filename = root+file(i);
    ret = importdata(filename);
    ROI(i) = {ret.ROI};
    bg_roi = ret.bgROI;
    
    %Access of max projection image matrix
    filename = root+string('modified_stack'+string(ids(i))+'.mat');
    load(filename);
    % imagesc(max(max_stack,[],3));
    img(i) = {max_stack};
    % hold on;
end

% Enhance the image using a kernel
% enhancedImage1 = img1;
% enhancedImage2 = img2;
% for i = size(img1,3)
%     enhancedImage1(:,:,i) = imfilter(img1(:,:,i), kernel);
%    
% end
% for i = size(img2,3)
%      enhancedImage2(:,:,i) = imfilter(img2(:,:,i), kernel);
% end

for j=1:size(img,2)
    img_curr = img{j};
    ROI_curr = ROI{j};
    features1 = [];
    points1 = [];
    for i=1:size(img_curr,3)
        
        [features,points] = getFeaturesROI(img_curr(:,:,i),ROI_curr,'SIFT');
    %     [features,points] = getFeatures(enhancedImage1(:,:,i),'SIFT');
        features1 = cat(1,features1,features);
        points1 = cat(1,points1,points);

    end
    feature_arr(j)={features1'};
    points_arr(j)={points1(:,1:2)'};
end

disp('extracted features');

disp('finding features');
%Cyclic consistency part


for i=1:size(img,2)
    
    views(i)=struct('desc',feature_arr{i},'img',max(img{i},[],3),'frame', ...
        points_arr{i},'nfeature',size(feature_arr{i},2),'filename',file{i});
end
pMatch = runGraphMatchBatch(views,'all',[],'wEdge', 0);
C = [];
for i = 1:length(views)
    cnt(:,i) = sum(views(i).frame,2)/double(views(i).nfeature);
    C = [C,views(i).frame - repmat(cnt(:,i),1,views(i).nfeature)];
end
%% Multi-Object Matching
% methods to try:
%  - 'pg': the proposed method, 
%         [Multi-Image Semantic Matching by Mining Consistent Features, CVPR 2018]
%  - 'spectral': Spectral method,
%          [Solving the multi-way matching problem by permutation synchronization, NIPS 2013]
%  - 'matchlift': MatchLift,
%          [Near-optimal joint object matching via convex relaxation, ICML 2014]   
%  - 'als': MatchALS,
%         [Multi-Image Matching via Fast Alternating Minimization, CVPR 2015]   
[jMatch,jmInfo] = runJointMatch(pMatch,C,'Method','als','univsize',10, ...
    'rank',3,'lambda', 1);
% save(savefile,'-append','jMatch','jmInfo');
%% Evaluate
X1 = pMatch2perm(pMatch); % pairwise matching result
X2 = pMatch2perm(jMatch); % joint matching result
n_img = length(imgList);
n_pts = length(X1)/n_img;
X0 = sparse(repmat(eye(ceil(n_pts)),n_img,n_img)); %groundtruth
% evaluate [overlap, precision, recall]
[o1,p1,r1] = evalMMatch(X1,X0);
[o2,p2,r2] = evalMMatch(X2,X0);
%% Visualize
if showmatch
    %view pairwise matches
    for i = 1:size(pMatch,1)
        for j = i+1:size(pMatch,2)
            clf;
            if ~isempty(pMatch(i,j).X)
                subplot('position',[0 0.5 1 0.48]);
                visPMatch(datapath,pMatch(i,j),3,'th',0.01);
                subplot('position',[0 0 1 0.48]);
                visPMatch(datapath,jMatch(i,j),3,'th',0.01);
                fprintf('%d-%d\n',i,j);
                pause
            end
        end
    end
end
% % matches = matchFeatures(features1,features2,'MatchThreshold',5.0,'Method','Approximate');
% [matches, scores] = vl_ubcmatch(features1', features2',1.3) ;
% matchedPoints1 = points1(matches(1,:),1:2);
% matchedPoints2 = points2(matches(2,:),1:2);
% I1 = max(img1,[],3);
% I2 = max(img2,[],3);
% thresh = 30;
% [matchedPoints1,matchedPoints2] = trimPoints(matchedPoints1,matchedPoints2,thresh);
% % figure; ax = axes;
% % showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2,'Parent',ax);
% % 
% 
% tform  = estimateGeometricTransform(matchedPoints1,matchedPoints2,'affine','MaxDistance',0.2,'Confidence',70);
% outputView = imref2d(size(I1));
% Ir = imwarp(I1,tform,'OutputView',outputView);
% % [H corrPtIdx] = findHomography(matchedPoints1',matchedPoints2');
% % tform = affine2d(round(H'));
% % [I21,I2_ref] = imwarp(I2,tform); % reproject img2
% % I21 = imtranslate(I2,-H(1:2,3)');
% time_elapsed = toc;
% figure;
% showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);
% figure;
% showMatchedFeatures(Ir,I2,matchedPoints1(1,:),matchedPoints2(1,:));
