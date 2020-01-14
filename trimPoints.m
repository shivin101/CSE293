function [points1,points2]= trimPoints(matchedPoints1,matchedPoints2,thresh)
   idx = abs(sum((matchedPoints1-matchedPoints2),2))<thresh
   points1 = matchedPoints1(idx,:)
   points2 = matchedPoints2(idx,:)
end