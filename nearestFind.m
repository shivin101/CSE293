function [matchNN,matchDist,vis] = nearestFind(iarea,dist,vis,idx,matchNN,matchDist)
    overlap = iarea(idx,:);
    euclid = dist(idx,:);
%     overlap(idx)=0;
%     euclid(idx)=9999;
    
    overlap = 1-(overlap/max(max(overlap),1));
    euclid = euclid/max(euclid);
    
    distVal = euclid+overlap;
    
    [B,I] = sort(distVal);
    [matchNN,matchDist,vis] = matchNeighbor(B,I,matchNN,matchDist,vis,idx);
    
end