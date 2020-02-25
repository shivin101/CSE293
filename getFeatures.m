function [features,points]= getFeatures(img,feature_type)

% img = linearize_val(img);
% img = linearize_val(normalize(double(img)));
% disp(max(img(:)));
% disp(min(img(:)));



sample_factor=200;
% disp(size(mid_p));
if feature_type=='SIFT'
    binSize = 25 ;
    magnif = 3 ;
%     Is = vl_imsmooth(single(img), sqrt((binSize/magnif)^2 - .25)) ;
    
    Is = zscore(single(img));
  
    [points,features] = vl_dsift(single(Is), 'size', binSize) ;
    points = downsample(points',sample_factor);
    features = downsample(features',sample_factor);
else
    [features,points]=extractFeatures(img,uint32(mid_p),'Method',feature);
end
end