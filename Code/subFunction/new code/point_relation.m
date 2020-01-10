function weight = point_relation(A,B)
global datas
[x1,y1,z1] = ind2sub(datas,A);
[x2,y2,z2] = ind2sub(datas,B);
% weight = 1
if abs(x1-x2)+abs(y1-y2)+abs(z1-z2)==1
    weight = 1;
elseif abs(x1-x2)+abs(y1-y2)+abs(z1-z2)==2
    weight = sqrt(2);
else
    weight = sqrt(3);
    %weight = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
end
    