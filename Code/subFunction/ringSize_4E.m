function s = ringSize_4E(ringdots)
A = ringdots(1,:);
B = ringdots(2,:);
C = ringdots(3,:);
D = ringdots(4,:);

a = vecdist(A,B);
b = vecdist(B,C);
c = vecdist(C,D);
d = vecdist(D,A);
e = vecdist(A,C);
s1 = sqrt((a+b+e)*(-a+b+e)*(a-b+e)*(a+b-e))/4;
s2 = sqrt((c+d+e)*(-c+d+e)*(c-d+e)*(c+d-e))/4;
s = s1+s2;
