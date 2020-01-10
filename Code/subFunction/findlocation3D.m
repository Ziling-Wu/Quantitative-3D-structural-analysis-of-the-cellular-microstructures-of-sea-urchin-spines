function loc = findlocation3D(x,y)
loc1 = find(x(1)==y(:,1));
loc2 = find(x(2)==y(:,2));
loc3 = find(x(3)==y(:,3));
loc4 = intersect(loc1,loc2);
loc = intersect(loc4,loc3);
