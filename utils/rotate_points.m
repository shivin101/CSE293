function rotated_points = rotate_points(points,y,x,theta)
% theta = theta/180;
center = [x,y];
R = [cosd(theta) sind(theta);-sind(theta) cosd(theta)]
% do the rotation...
s = points - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s';           % apply the rotation about the origin
vo = so' + center;   % shift again so the origin goes back to the desired center of rotation
% this can be done in one line as:
% vo = R*(v - center) + center
% pick out the vectors of rotated x- and y-data
rotated_points = vo;