clearvars

newskel = zeros(100,100,100);
newskel(1:20:100,1:20:100,1:25:100)=1;
%%
datas = [100,100,100];
global datas
cropSkel_node = find(newskel);
% D = graph(cropSkel_node);
%%
noise = -5:5;
[x,y,z] = ind2sub(datas,cropSkel_node);
for i = 1:length(x)
    if (x(i)>10)&&(x(i)<90)
        newx = x(i)+noise(randperm(length(noise),1));
        x(i) = newx;
    elseif (y(i)>10)&&(y(i)<90)
        newy = y(i)+noise(randperm(length(noise),1));
        y(i) = newy;
    elseif (z(i)>10)&&(z(i)<90)
        newz = z(i)+noise(randperm(length(noise),1));
        z(i) = newz;
    end
end
newnode = sub2ind(datas,x,y,z);
PlotScatter(newnode(:,1),20,'r');
hold on
PlotScatter(cropSkel_node(:,1),20,'g');
FigureFormat
cropSkel_node = newnode;
%%
clear iconnect
[x,y,z] = ind2sub(datas,cropSkel_node);
xyz = [x,y,z];
out = squareform(pdist(xyz));
% iconnect = zeros(1,2);
iconnect(:,1) = repelem(1:length(cropSkel_node),5);
% connect each 
for i = 1:length(cropSkel_node)
    loc = zeros(6,1);
    dis = zeros(6,1);
    val = out(i,:);
    for j=1:6
        [flag,idx] = min(val);
        loc(j) = idx;
        dis(j) = flag;
      % remove for the next iteration the last smallest value:
        val(idx) = 100;
    end
    iconnect((iconnect(:,1)==i),2) = (loc(2:6));
    iconnect((iconnect(:,1)==i),3) = (dis(2:6));
end
%% angle
for j = 1:length(iconnect)
    [x1,y1,z1] = ind2sub([100,100,100],cropSkel_node(iconnect(j,1)));
    [x2,y2,z2] = ind2sub([100,100,100],cropSkel_node(iconnect(j,2)));
    iconnect(j,4)=atand(sqrt((y1-y2).^2 + (z1-z2).^2)/abs(x1-x2));
end
%%
for i = 1:length(cropSkel_node)
    loc1 = (find(iconnect(:,1)==i));%points index connect to i
    loc1_ind = iconnect(loc1,2);%points index connect to i
    loc2 = find(iconnect(:,2)==i);%points index connect to i
    loc2_ind = iconnect(loc2,1);%points index connect to i
    loc = intersect(loc1_ind,loc2_ind)
    if length(loc)>=3
        c1 = setdiff(loc1_ind,loc);
        d1 = [];
        if ~isempty(c1)
           for j = 1:length(c1) 
                d1(j) = loc1(find(loc1_ind==c1(j)));
           end
        end
        c2 = setdiff(loc2_ind,loc);
        d2 = [];
        if ~isempty(c2)
            
           for j = 1:length(c2) 
                d2(j) = loc2(find(loc2_ind==c2(j)));
           end
        end
        iconnect(d1,:)=[];
        iconnect(d2,:)=[];
    end
        
end
%%
iconnect_sort = sort(iconnect(:,1:2),2);
iconnect_sort = [iconnect_sort,iconnect(:,3:end)];
iconnect_sort = unique(iconnect_sort,'rows');
%% no 3-edge ring
for i = 1:size(cropSkel_node,1)
    loc1 = (find(iconnect_sort(:,1)==i));%points index connect to i
    loc1_ind = iconnect_sort(loc1,2);%points index connect to i
    loc2 = find(iconnect_sort(:,2)==i);%points index connect to i
    loc2_ind = iconnect_sort(loc2,1);%points index connect to i
    loc_ind = [loc2_ind;loc1_ind];
    loc = [loc2;loc1];
    for j = 1:length(loc)-1
        for k = 2:length(loc)
            ring3 = find((iconnect_sort(:,1)==min(loc_ind(j),loc_ind(k)))&(iconnect_sort(:,2)==max(loc_ind(j),loc_ind(k))));
            if ~isempty(ring3)
                ind = [min(loc(j),loc(k)),max(loc(j),loc(k)),ring3];
                dis = [iconnect_sort(loc(j),3),iconnect_sort(loc(k),3),iconnect_sort(ring3,3)];
                iconnect_sort(ind(dis == max(dis)),:)=[];
            end
%             dis == max(dis)
        end
        
    end
end
%%
n_node = size(cropSkel_node,1);
idx=randperm(n_node);
n3 = sort(idx(1:round(n_node*0.7)));
n4 = sort(idx(round(n_node*0.7)+1:round(n_node*0.95)));
n5 = sort(idx(round(n_node*0.95)+1:end));
cropSkel_node(n3,2)=3;
cropSkel_node(n4,2)=4;
cropSkel_node(n5,2)=5;
for i = 1:n_node
    loc1 = (find(iconnect_sort(:,1)==i));%points index connect to i
    loc1_ind = iconnect_sort(loc1,2);%points index connect to i
    loc2 = find(iconnect_sort(:,2)==i);%points index connect to i
    loc2_ind = iconnect_sort(loc2,1);%points index connect to i
    loc_ind = [loc2_ind;loc1_ind];
    loc = [loc2;loc1];
    if length(loc)>cropSkel_node(i,2)
        criterion = iconnect_sort(loc,4);
        toBeDelete = length(loc)-cropSkel_node(i,2)
        ind = sortrows([criterion,loc],'descend')
        iconnect_sort(ind(1:toBeDelete,2),:)=[];
    end
end
%%
iconnect_sort(:,1) = cropSkel_node(iconnect_sort(:,1));  
iconnect_sort(:,2) = cropSkel_node(iconnect_sort(:,2));  
%%

clf
[x,y,z] = ind2sub(datas,iconnect_sort(:,1));
xyz1 = [x,y,z];
[x,y,z] = ind2sub(datas,iconnect_sort(:,2));
xyz2 = [x,y,z];
for i = 1:size(xyz1,1)
    xx = [xyz1(i,1),xyz2(i,1)];
    yy= [xyz1(i,2),xyz2(i,2)];
    zz= [xyz1(i,3),xyz2(i,3)];
    plot3(yy,xx,zz,'b','LineWidth',1)
%     hold on
%     plot3(,'r')
    hold on
end
PlotScatter(xyz1)
FigureFormat
%%
idx=randperm(length(x));
n3 = sort(idx(1:round(length(x)*0.5)));
n4 = sort(idx(round(length(x)*0.5)+1:round(length(x)*0.9)));
n5 = sort(idx(round(length(x)*0.9)+1:end));
%%
criteria = 100;
%%
