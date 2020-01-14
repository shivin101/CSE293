X = rand(2,50)*500 - [0;50];
Y = [cosd(5) sind(5);-sind(5) cosd(5)]*(X-[300;0])+[300;0];
figure;scatter(X(1,:),X(2,:))
hold on
scatter(Y(1,:),Y(2,:))
for i=1:50
    line([X(1,i) Y(1,i)] , [X(2,i), Y(2,i)],'Color','k');
end
axis ij
axis equal
%%
X = rand(2,100)*300 - [0;150];
Y = [cosd(2) sind(2);-sind(2) cosd(2)]*(X-[1000;0])+[1000;0];
figure;scatter(X(1,:),X(2,:))
hold on
scatter(Y(1,:),Y(2,:))
for i=1:100
    line([X(1,i) Y(1,i)] , [X(2,i), Y(2,i)],'Color','k');
end