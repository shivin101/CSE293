clc;
clear all;
close all;

tic
startup;
ids=[1,2];

root = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\New_data\';
dir_list = dir(root);
dir_list = dir_list(3:end);

dir_idx = 10;
val = dir_list(dir_idx);
mod_dir = [val.folder,'\',val.name];
ROI ={1,size(ids,2)};

num_sessions = dir([mod_dir,'\','*.mat']);

file_list = {};
%Hard coding some naming conventions in this stage of the process
for i = 1:size(ids,2)
    file_list(i) = {[mod_dir,'\','modified_stack',num2str(ids(i)),'.mat']};
    load(file_list{i});
    
    %Create an array of ROI from preprocessed arrays
%     filename = root+file(i);
%     ret = importdata(filename);
    Roi(i) = {ROI};
%     bg_roi = ret.bgROI;
    
    %Create an array of Image array
    
    % imagesc(max(max_stack,[],3));
    img(i) = {max_stack};
    % hold on;
end

%ROI is list of cells ie 1xN where N is the number of sessions
%Each ROI is a cell of matrices

%img is a cell of images 
%Each image is already preprocessed with size of the format HxWxD


%Extract features based on ROI centers
for j=1:size(img,2)
    img_curr = img{j};
    ROI_curr = Roi{j};
    
    %Just place holders for this part of the code
    features1 = [];
    points1 = [];
    for i=1:size(img_curr,3)
        
        [features,points] = getFeaturesROI(img_curr(:,:,i),ROI_curr,'SIFT');
        features1 = cat(1,features1,features);
        points1 = cat(1,points1,points);

    end
    feature_arr(j)={features1};
    points_arr(j)={points1(:,1:2)};
end


disp('extracted features');
%Just a placeholder for pairs of images
features1 = feature_arr{1};
features2 = feature_arr{2};
points1 = points_arr{1};
points2 = points_arr{2};
img1 = img{1};
img2 = img{2};

%Find matches between descriptors
[matches, scores] = vl_ubcmatch(features1', features2',1.1) ;
matchedPoints1 = points1(matches(1,:),1:2);
matchedPoints2 = points2(matches(2,:),1:2);

%Create a max projection of the entire image
I1 = max(img1,[],3);
I2 = max(img2,[],3);
thresh = 30;

%Trim points from matches points based on a threshold
[matchedPoints1,matchedPoints2] = trimPoints(matchedPoints1,matchedPoints2,thresh);

%Use RANSAC to find a geometric transformation
tform  = estimateGeometricTransform(matchedPoints1,matchedPoints2,...
    'affine','MaxDistance',0.2,'Confidence',90);

%Apply transformation
outputView = imref2d(size(I1));
Ir = imwarp(I1,tform,'OutputView',outputView);
time_elapsed = toc;
figure;
showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);
figure;
showMatchedFeatures(Ir,I2,matchedPoints1(1,:),matchedPoints2(1,:));
