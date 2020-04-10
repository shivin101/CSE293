function [pair_ROI]=findPairROI(tform,ROI1,ROI2,midp1,midp2)
    %%Transform the ROI coordinates to the new coordinate frame
    % using the transformation obtained from RANSAC
    for i=1:size(ROI1,2)
        points=ROI1{i};
        u=points(:,1);
        v=points(:,2);
        [x,y]=transformPointsForward(tform,u,v);
        transform_points = cat(2,x,y);
        ROI1{i}=transform_points;
    end
    
    %%Calcualate the distance metrics for your ROI's 
    %In this case we use the geometric distance between 
    IOU = zeros(size(ROI1,2),size(ROI2,2));
    dist = pdist2(midp1,midp2);
    iarea = zeros(size(ROI1,2),size(ROI2,2));
    for i = 1:size(ROI1,2)
        for j = 1:size(ROI2,2)
%             ihull = intersectionHull('vert',ROI1{i},'vert',ROI2{j});
            ipoly = intersect(polyshape(ROI1{i}),polyshape(ROI2{j}));
            iarea(i,j) = area(ipoly);
%             iarea(j,i) = iarea(i,j);
        end
    end
    vis = zeros(size(ROI2,2),1);
    matchNN=zeros(size(ROI1,2),1);
    matchDist=ones(size(ROI1,2),1)*9999;
    for idx=1:size(ROI1,2)
        [matchNN,matchDist,vis]=nearestFind(iarea,dist,vis,idx,matchNN,matchDist);
    end
    pair_ROI = matchNN;
;end