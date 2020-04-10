clc;
clear all;
close all;

tic

ids=[1,2,3,4];
denoise = 0;
root = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\New_data\';
dir_list = dir(root);
dir_list = dir_list(3:end);

dir_idx =10;
val = dir_list(dir_idx);
mod_dir = [val.folder,'\',val.name];
%ROI ={1,size(ids,2)};

num_sessions = dir([mod_dir,'\','*.mat']);

file_list = {};
%Hard coding some naming conventions in this stage of the process
for i = 1:size(ids,2)
    file_list(i) = {[mod_dir,'\','modified_stack',num2str(ids(i)),'.mat']};
    file(i) = file_list(i);
    load(file_list{i});
    
    %Create an array of ROI from preprocessed arrays
%     filename = root+file(i);
%     ret = importdata(filename);
    Roi(i) = {ROI};
%     bg_roi = ret.bgROI;
    
    %Create an array of Image array
    
    % imagesc(max(max_stack,[],3));
    img(i) = {max_stack};
    if denoise
        for j =1:size(img{i},3)
            curr_img = max_stack;
            den_img = preprocess_data(double(curr_img(:,:,j)));
            curr_img(:,:,j)=den_img;
        end
        img(i)={curr_img};
    end
    % hold on;
end


%*********Dead Code***********
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


%Extract features based on ROI centers
for j=1:size(img,2)
    img_curr = img{j};
    ROI_curr = Roi{j};
    features1 = [];
    points1 = [];
    for i=1:size(img_curr,3)
        
        [features,points,mid_points] = getFeaturesROI(img_curr(:,:,i),ROI_curr,'SIFT');
    %     [features,points] = getFeatures(enhancedImage1(:,:,i),'SIFT');
        features1 = cat(1,features1,features);
        points1 = cat(1,points1,points);

    end
    feature_arr(j)={features1'};
    midp_arr(j)={mid_points(:,1:2)};
    points_arr(j)={points1(:,1:2)'};
end


disp('extracted features');

disp('finding features');


%Cyclic consistency part
for i=1:size(img,2)
    
    views(i)=struct('desc',feature_arr{i},'img',uint8(max(img{i},[],3)/3),'frame', ...
        points_arr{i},'nfeature',size(feature_arr{i},2));
end
pMatch = runGraphMatchBatch(views,file,'all',[],'wEdge', 0);
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
[jMatch,jmInfo] = runJointMatch(pMatch,C,'Method','pg','univsize',20, ...
    'rank',3,'lambda', 1);
% save(savefile,'-append','jMatch','jmInfo');
%% Evaluate
% X1 = pMatch2perm(pMatch); % pairwise matching result
% X2 = pMatch2perm(jMatch); % joint matching result
% n_img = length(img);
% n_pts = length(X1)/n_img;
% X0 = sparse(repmat(eye(ceil(n_pts)),n_img,n_img)); %groundtruth
% % evaluate [overlap, precision, recall]
% [o1,p1,r1] = evalMMatch(X1,X0);
% [o2,p2,r2] = evalMMatch(X2,X0);
showmatch = true;
%% Visualize
if showmatch
    %view pairwise matches
    for i = 1:size(pMatch,1)
        for j = i+1:size(pMatch,2)
            clf;
            if ~isempty(pMatch(i,j).X)
                subplot('position',[0 0.5 1 0.48]);
                visPMatch(pMatch(i,j),views([i,j]),1,'th',0.01);
                subplot('position',[0 0 1 0.48]);
                visPMatch(jMatch(i,j),views([i,j]),1,'th',1.0);
                fprintf('%d-%d\n',i,j);
                pause
            end
        end
    end
end
roiDict = cell(size(pMatch,1),size(pMatch,2))
 for i = 1:size(pMatch,1)
        for j = i+1:size(pMatch,2)
            tformArr(i,j)={getTform(jMatch(i,j),views([i,j]))};
            fprintf('%d-%d\n',i,j);
            
 
           roiDict(i,j) = {findPairROI(tformArr{i,j},Roi{i},Roi{j},midp_arr{i},midp_arr{j})};
        end
 end
