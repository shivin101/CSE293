clc;
clear all;
close all;

id1 = 3;
id2 = 6;
root = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\M1-2\';
file = 'KA0013_1801'+string(num2str(id1+5,'%02d'))+'_001_001_summed_50.roi';

filename = root+file;
ret = importdata(filename);
ROI1 = ret.ROI;
bg_roi = ret.bgROI;

filename = root+string('modified_stack'+string(id1)+'.mat');
load(filename);
% imagesc(max(max_stack,[],3));
img1 = max_stack;
% hold on;

root = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\M1-2\';
file = 'KA0013_1801'+string(num2str(id2+5,'%02d'))+'_001_001_summed_50.roi';

filename = root+file;
ret = importdata(filename);
ROI2 = ret.ROI;
bg_roi = ret.bgROI;

filename = root+string('modified_stack'+string(id2)+'.mat');
load(filename);
% imagesc(max(max_stack,[],3));
img2 = max_stack;


    
[points1,features1] = maximal_features(ROI1,img1);

features2 = [];
points2 = [];

[points2,features2] = maximal_features(ROI2,img2);





% matches = matchFeatures(features1,features2,'MatchThreshold',5.0,'Method','Approximate');
[matches, scores] = vl_ubcmatch(features1', features2',1.3) ;
matchedPoints1 = points1(matches(1,:),1:2);
matchedPoints2 = points2(matches(2,:),1:2);
I1 = max(img1,[],3);
I2 = max(img2,[],3);
thresh = 30;
[matchedPoints1,matchedPoints2] = trimPoints(matchedPoints1,matchedPoints2,thresh);
figure; ax = axes;
showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2,'Parent',ax);
% for i=1:size(ROI,2)
%    points=ROI{i};
%    x= points(:,1);
%    y= points(:,2);
%    k=convhull(x,y);
%    plot(x(k),y(k),'r-');
% end
% 
% mid_p = [];
% for i=1:size(ROI,2)
%    points=ROI{i};
%    mid_p = [mid_p;mean(points)];
%    
% end
% plot(mid_p(:,1),mid_p(:,2),'g+');
% 
% 
% extractFeatures(img,uint32(mid_p),'Method','SURF');
