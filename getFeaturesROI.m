function [features,points,mid_p]= getFeaturesROI(img,ROI,feature_type)

% img = linearize_val(img);
% img = linearize_val(normalize(double(img)));
% disp(max(img(:)));
% disp(min(img(:)));

mid_p = [];
feat_size = 11;
for i=1:size(ROI,2)
   points=ROI{i};
   mid_p = [mid_p;[mean(points),feat_size,0]];
   
end
% disp(size(mid_p));
if feature_type=='SIFT'
    [points,features]=vl_sift(single(img),'frames',mid_p','orientations');
    points = points';
    features = features';
else
    [features,points]=extractFeatures(img,uint32(mid_p),'Method',feature);
end
end