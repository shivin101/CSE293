function [matchedPoints1,matchedPoints2] = get_corresp(img1,ROI1,img2,ROI2)
features1 = [];
points1 = [];
for i=1:size(img1,3)
    
    [features,points] = getFeaturesROI(img1(:,:,i),ROI1,'SIFT');
%     [features,points] = getFeatures(img1(:,:,i),'SIFT');
    features1 = cat(1,features1,features);
    points1 = cat(1,points1,points);

end
disp('extracted features 1');
features2 = [];
points2 = [];

for i=1:size(img2,3)
    
    [features,points] = getFeaturesROI(img2(:,:,i),ROI2,'SIFT');
%     [features,points] = getFeatures(img2(:,:,i),'SIFT');
    features2 = cat(1,features2,features);
    points2 = cat(1,points2,points);

end
disp('extracted features 2');

disp('finding corresp features');
disp(size(features1));
% matches = matchFeatures(features1,features2,'MatchThreshold',5.0,'Method','Approximate');
[matches, scores] = vl_ubcmatch(features1', features2',1.1) ;
matchedPoints1 = points1(matches(1,:),1:2);
matchedPoints2 = points2(matches(2,:),1:2);
I1 = max(img1,[],3);
I2 = max(img2,[],3);
thresh = 30;
[matchedPoints1,matchedPoints2] = trimPoints(matchedPoints1,matchedPoints2,thresh);
end