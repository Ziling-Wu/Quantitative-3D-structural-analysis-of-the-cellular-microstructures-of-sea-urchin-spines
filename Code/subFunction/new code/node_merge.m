function [node_new,bran_new] = node_merge(skel)
%merge nodes and recombine the branches 

%% save branches into cell
clear bran
length_threshold = 2; % branch length threshold
xsize = size(skel,1);ysize = size(skel,2);zsize = size(skel,3);
index =1;
new_BW = zeros(xsize,ysize,zsize);
for iter = 1:max(BW(:))
    corbran = find(BW == iter);
    [cora,corb,corc] = ind2sub([xsize ysize zsize],corbran); % find nonzero component
    element = [cora,corb,corc];
    if size(element,1) >= length_threshold
       bran{index,1} = element;
       index = index+1;
    end
end
%% save branches into cell-calculate thickness
clear new_bran
length_threshold = 3; % branch length threshold
xsize = size(skel,1);ysize = size(skel,2);zsize = size(skel,3);
index =1;
new_BW = zeros(xsize,ysize,zsize);
for iter = 1:max(BW(:))
    corbran = find(BW == iter);
    [cora,corb,corc] = ind2sub([xsize ysize zsize],corbran); % find nonzero component
    element = [cora,corb,corc];
    if size(element,1) >= length_threshold
       new_bran{index,1} = element;
       index = index+1;
    end
end
%%
for i =1:size(bran,1)
    clear element
    element = bran{i};
    for j = 1:size(element,1)
        new_BW(element(j,1),element(j,2),element(j,3))=BW(element(j,1),element(j,2),element(j,3));
    end
end
%% connected nodes to these branches
numbran = size(new_bran,1);
 
%check which nodes belongs to these branches, if there is no node, just put
%into the endpoint

connectedbran = zeros(numbran,6);
for i = 1:numbran
    index=1;% one end
    element = new_bran{i};
    % check every element's nearby points to see if it node
    numele = size(element,1);
    for k = 1:numele
            [val,s] = survalue(element(k,:),skel,xsize,ysize,zsize);
            loc = find(val);
            numin = length(loc); % numin=surrounding branches+1
            if numin > 2 % it has a neighbor
            for j = 1:numin
                cora = s(loc(j),1);corb = s(loc(j),2);corc = s(loc(j),3);
                cornode = [cora,corb,corc];
                flag = ismember(cornode, node,'rows');
                if flag == 1
                    connectedbran(i,index:index+2) = cornode;
                    index = index+3;
                end
            end
            elseif numin == 2 % endpoint
                cora = element(k,1);corb = element(k,2);corc = element(k,3);
                cornode = [cora,corb,corc];
                connectedbran(i,index:index+2) = cornode;
                index = index+3;
            end
    end 
end
%% check start point and endpoint for each branch
% here only consider branch without node
% start and endpoint will only have one neighbor
numbran = length(new_bran);
startbran = zeros(numbran,6);

% check start point and endpoint
for i = 1:numbran
    branmatrix = zeros(xsize,ysize,zsize); % construct matrix only contain one branch
    index = 1;
    element = new_bran{i};
    numele = size(element,1);
    for t =1:numele
        branmatrix(element(t,1),element(t,2),element(t,3))=1;
    end
    for k = 1:numele
            [val,s] = survalue(element(k,:),branmatrix,xsize,ysize,zsize);
            loc = find(val);
            numin = length(loc); % surrounding branches+1
            if numin == 2 %start or end point
                cora = element(k,1);corb = element(k,2);corc = element(k,3);
                cornode = [cora,corb,corc];
                startbran(i,index:index+2)=cornode;
                index = index+3;
            end
    end
end


%% some nodes may not coonected to a branch
clear new_node
index = 1;
P1 = connectedbran(:,1:3);
P2 = connectedbran(:,4:6);
for iter = 1 : size(node,1)
    replace1 = find(node(iter,1) == P1(:,1) & node(iter,2) == P1(:,2) & node(iter,3) == P1(:,3));
    replace2 = find(node(iter,1) == P2(:,1) & node(iter,2) == P2(:,2) & node(iter,3) == P2(:,3));
    replace = [replace1; replace2];
    if ~isempty(replace)
        new_node(index,1:3) = node(iter,:);
        index = index+1;
    end
end
node_new = new_node;
%% 
numnode = size(new_node,1);
index_out = 1;
index_in = 1;
clear node_out
node_in = zeros(1,3);
nodetype = cell(1,1);
threhold_merge = 4;
index = 1;
%node_in(1,:)=node_new(1,:);
for iter = 1:numnode-1
    loc = findlocation3D(node_new(iter,:),node_in);%
    if isempty(loc)
       freq = branch_connected(node_new(iter,1:3),new_BW,skel);
       array = iter+1:numnode;
       distance = merge_thres;
       clear temp
       for i = 1:length(array)/2
           flag = abs(node_new(array(i),1)-node_new(iter,1))<threhold_merge&abs(node_new(array(i),2)-node_new(iter,2))<threhold_merge&...
               abs(node_new(array(i),3)-node_new(iter,3))<threhold_merge;
           if flag
               temp_distance = vecdist(node_new(array(i),:),node_new(iter,1:3));
               if temp_distance<distance
                   distance = temp_distance;
                   temp = node_new(array(i),:);
                   if distance == 1
                       break
                   end
               end
           end
       end
       if distance ~= merge_thres
           loc1 = findlocation3D(temp,node_in);
           if isempty(loc1)
               node_in(index,1:3) = temp;
               freq1 = branch_connected(temp,new_BW,skel);
               freq = [freq,freq1];
               nodetype{index} = freq;
               index = index+1;
           else
               freq1 = nodetype{loc1};
               freq = [freq,freq1];
               nodetype{loc1} = freq;
           end
       else%don't find neighbor
           node_in(index,1:3) = node_new(iter,1:3);
           nodetype{index} = freq;
           index = index+1;
       end
    else % ~isempty(loc)
        freq = nodetype{loc};
        array = iter+1:numnode;
        distance = merge_thres;
        clear temp
       for i = 1:length(array)
           flag = abs(node_new(array(i),1)-node_new(iter,1))<threhold_merge&abs(node_new(array(i),2)-node_new(iter,2))<threhold_merge&...
               abs(node_new(array(i),3)-node_new(iter,3))<threhold_merge;
           if flag
               temp_distance = vecdist(node_new(array(i),:),node_new(iter,1:3));
               if temp_distance<distance
                   distance = temp_distance;
                   temp = node_new(array(i),:);
                   if distance == 1
                       break
                   end
               end
           end
       end
       if distance ~= merge_thres
           loc1 = findlocation3D(temp,node_in);   
           if isempty(loc1)
               node_in(loc,1:3) = temp;
               freq1 = branch_connected(temp,new_BW,skel);
               freq = [freq,freq1];
               nodetype{loc} = freq;
           else
               freq1 = nodetype{loc1};
               freq = [freq,freq1];
               nodetype{loc1} = freq;
               nodetype{loc}=[];
               node_in(loc,1:3)=0;
           end
           
       else 
           continue
       end
    end
end
%%
output = node_in;
clear node_in
index=1;
for j = 1:size(nodetype,2)
    output(j,4) = length(nodetype{j})-2*(length(nodetype{j})-length(unique(nodetype{j})));
    if output(j,4) <3
        output(j,4)=0;
    end
end
for j = 1:size(nodetype,2)
    if output(j,4) ~= 0
        node_in(index,:) = output(j,:);
        output_orientation{index} = nodetype{j};
        index=index+1;
    end
end
