function [connectTemp,newnode,branTemp] = merge_Allnear(connectTemp,nodeTemp,branTemp,dist)
% dist = 10;
% connectTemp=c2;
% branTemp=b2;
% nodeTemp = n2;
global datas
%% Calculate distance for every pair 
[x1,y1,z1] = ind2sub(datas,nodeTemp(:,1));
X = [x1,y1,z1];
D = pdist2(X,X);
%%
nodeOutId =[];
for i = 1:size(nodeTemp,1)
    flag = ismember(i,nodeOutId);
    if flag
        continue
    else
        nodeClusterId = i;
        ind = find((D(i,:)<=dist)&(D(i,:)>0));%find nearby points
        if ~isempty(ind) % exist nearby points
            nodeClusterId = [nodeClusterId,ind];
            NearNodeId = findNearNode(ind,D,nodeClusterId,dist);
            while ~isempty(NearNodeId)
                nodeClusterId = [nodeClusterId,NearNodeId];
                NearNodeId = findNearNode(NearNodeId,D,nodeClusterId,dist);
            end
        end
        if length(nodeClusterId)>1
            nodeCluster = nodeTemp(nodeClusterId,1);
            [x,y,z] = ind2sub(datas,nodeCluster);
            center = sub2ind(datas,round(mean(x)),round(mean(y)),round(mean(z)));
            for j = 1:length(nodeClusterId)
                cc = find(connectTemp(:,1) == nodeCluster(j));
                cd = find(connectTemp(:,2) == nodeCluster(j));
                if ~isempty(cc)
                   connectTemp(cc,1) = center;
                end
                if ~isempty(cd)
                   connectTemp(cd,2) = center;
                end
            end
        end
        nodeOutId = [nodeOutId,nodeClusterId];
    end
end
 
%%
cluster = find(connectTemp(:,1)==connectTemp(:,2));
cluster = flip(cluster);
for k = 1:length(cluster)
    newconnect = [connectTemp(1:cluster(k)-1,:);connectTemp(cluster(k)+1:end,:)];
    clear connectTemp; 
    connectTemp = newconnect;
    branTemp{cluster(k)}=[];
end
id = cellfun('length',branTemp);branTemp(id==0)=[];
newnode = unique(connectTemp);
for i = 1:length(newnode)
    numcon = length(find(connectTemp(:,1)==newnode(i,1)))+length(find(connectTemp(:,2)==newnode(i,1)));
    newnode(i,2) = numcon;
end

        