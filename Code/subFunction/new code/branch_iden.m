function [bran,connect,node]= branch_iden(skel,node)

global datas
%     node = node_identification(skel);
%     save('node.mat','node');
    nearnode = node(:,1);
    skel(ind2sub(datas,node(:,1)))=0;

BW = bwlabeln(skel,26);

indexbran = 1;
for iter = 1:max(BW(:))
    %sort the voxels for each branch
    corbran = find(BW == iter);
%     fprintf('Iteration %d how long the branch %d\n',iter,size(corbran));
    if size(corbran,1)>1  %for more than 1 voxel branches
        
         indexend=1;
         clear startnode
        for k = 1:size(corbran,1) %for all voxels on this branch
                n=neighbor(corbran(k),datas);
                numin=sum(double(ismember(n,corbran)));
                if numin == 2 % end point identification
                    temp=(double(ismember(n,nearnode)));
                        if sum(temp)==1 % connected to one and only one node
                            startnode(indexend)=n(temp==1); % register the node 
                            indexend=indexend+1;
                        end
                end
        end
%         fprintf('connected node %d\n',startnode);

        if indexend==2+1 % both ends are connected to one and only one node
%             fprintf('%d\n',startnode);
            connect(indexbran,:)=startnode';% regsiter the nodes
            bran{indexbran} = corbran; % and the branch
            indexbran=indexbran+1;
            
        elseif indexend==1+1  % one ended branch
            
            %bran{indexbran}=[[startnode,0]';corbran];   % regsiter the one ended nodes and the branch
            %indexbran=indexbran+1;

        end
        
    else % for 1 voxel branches 
        
       indexend=1;clear startnode
        n=neighbor(corbran(1),datas);
                numin=sum(double(ismember(n,corbran))); 
                    temp=(double(ismember(n,nearnode)));
                        if sum(temp)==2 % this single voxel branch is connected to two nodes
                            startnode(indexend:indexend+1)=n(temp==1); % register the node 
                            indexend=indexend+1;
                        end
                        
        if indexend==1+2 % this single voxel branch is connected to two nodes
            
            connect(indexbran,:)=startnode';% regsiter the nodes
            bran{indexbran} = corbran; % and the branch
            indexbran=indexbran+1;

        end
%         fprintf('connected node %d\n',startnode);
    end
    
end
        
        

end