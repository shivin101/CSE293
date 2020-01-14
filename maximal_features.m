function [points,features] = maximal_features(ROI,img)
    selected_id = zeros(size(ROI,2),1);
    mid_p=[];
    selected_id=[];
    feat_size=11;
    
%   Select the point of maximal intensity in each of the volumes
    for i=1:size(ROI,2)
       points=ROI{i};
       point = uint32(mean(points));
       mid_p = [mid_p;[mean(points),feat_size,0]];
       [intensity,idx] = max(img(point(1),point(2),:));
       selected_id = [selected_id;idx];
    end
    
    
%   Extract the features from each of the slices
    points = [];
    features = [];
    
    for i=1:size(img,3)
        ids = (selected_id==i);
        selection = mid_p(ids,:,:,:);
        [point,feature]=vl_sift(single(img(:,:,i)),'frames',selection','orientations');
        points = [points;point'];
        features = [features;feature'];
    end
    
end