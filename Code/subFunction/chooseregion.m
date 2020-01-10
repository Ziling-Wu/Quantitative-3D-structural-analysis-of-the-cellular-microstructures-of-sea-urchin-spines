function loc = chooseregion(A,B,dis)
% A is nx3; B is 1x3 dis is the wanted distance between each element in B and A.
global datas
[x0, y0, z0] = ind2sub(datas,B);
loc1 = find(abs(A(:,1)-x0)<dis);
loc2 = find(abs(A(:,2)-y0)<dis);
loc3 = find(abs(A(:,3)-z0)<dis);
loc_t = intersect(loc1,loc2);
loc = intersect(loc_t,loc3);

