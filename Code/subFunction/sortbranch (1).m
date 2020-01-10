function [bran,node]= sortbranch(node, skel)

global datas
    skel(ind2sub(datas,node(:,1)))=0;

BW = bwlabeln(skel,26);

indexbran=1;
for iter = 1:max(BW(:))
    
    %sort the voxels for each branch
    corbran = find(BW == iter);
   
    % identify the starting and ending points and its connected nodes
     indexend=1;
    for k = 1:size(corbran,1) %for all voxels on this branch
            n=neighbor(corbran(k),datas);
            numin=sum(double(ismember(n,corbran)));
            if numin == 2 %start or end point
                temp=(double(ismember(n,node(:,1))));
                    if sum(temp)==1 %connected to one and only one node
                        startnode(indexend)=n(temp==1); % register the node 
                        indexend=indexend+1;
                    end
            end
    end

    if indexend==3 % two ends both connected to one and only one node
        
        bran{indexbran} = [startnode';corbran]; % regsiter the node and the branch
        indexbran=indexbran+1;
   
    end
    
end

    

end