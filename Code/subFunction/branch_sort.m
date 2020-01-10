function [bran]= branch_sort(node, skel)

global datas
skel(ind2sub(datas,node(:,1)))=0;
BW = bwlabeln(skel,26);
indexbran=1;

for iter = 1:max(BW(:))
    %sort the voxels for each branch
    corbran = find(BW == iter);
    if size(corbran,1)>1  %for more than 1 voxel branches
        indexend=0;
        for k = 1:size(corbran,1) %for all voxels on this branch
                n = neighbor(corbran(k),datas);
                numin = sum(double(ismember(n,corbran)));%
                if numin == 2 % end point identification
                    temp=(double(ismember(n,node(:,1))));
                    if sum(temp)==1 % connected to one and only one node
                        indexend=indexend+1;
                        startnode(indexend)=n(temp==1); % register the node
                    end
                end
        end
        if indexend==2 % both ends are connected to one and only one node
            bran{indexbran} = [startnode';corbran]; % regsiter the node and the branch
            indexbran=indexbran+1;
        end
    else % for 1 voxel branches
        indexend = 0;
        n = neighbor(corbran(1),datas);
        temp = (double(ismember(n,node(:,1))));
        if sum(temp) == 2 % this single voxel branch is connected to two nodes
            startnode=n(temp==1); % register the node 
            indexend=indexend+1;
        end
        if indexend==1 % this single voxel branch is connected to two nodes
            bran{indexbran} = [startnode';corbran]; % regsiter the node and the branch
            indexbran=indexbran+1;
        end
    end  
end
end