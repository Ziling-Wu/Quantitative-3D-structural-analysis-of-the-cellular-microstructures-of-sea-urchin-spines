function node = node_analysis(skel)
%get node positions
xsize = size(skel,1);ysize = size(skel,2);zsize = size(skel,3);

skelcor = find(skel);
[skelcorx,skelcory,skelcorz] = ind2sub([xsize ysize zsize],skelcor);

skelxyz = [skelcorx,skelcory,skelcorz];


lens =  size(skelxyz,1); %number of nozeros
index = 1;
for iter = 1:lens
c = skelxyz(iter,:);
[val,s] = survalue(c,skel,xsize,ysize,zsize);
nonzero = length(find(val));
if nonzero >= 4 % which should be a node
        node(index,1) = c(1);
        node(index,2) = c(2);
        node(index,3) = c(3);
        index = index+1;
end
end