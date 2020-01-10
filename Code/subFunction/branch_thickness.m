function [thickness,newbran_thickness,newconnect_thickness] = branch_thickness(bran,connect,node,new_img)
global datas
%%
% cont = zeros(datas);
% for i = 1:size(new_img,3)
%     cont(:,:,i) = edge(new_img(:,:,i));
% end
cont = edge3(new_img,'approxcanny',0.6);
%% newbran , newconnect middle point and direction
% establish new bran and new connect matrix since we only need to consider
% two-node connected branch. 
newbran = bran;
newconnect = connect;
% for i = 1:length(bran)
%     flag = ismember(connect(i,1),node)&ismember(connect(i,2),node);
%     if flag
%         newconnect = [newconnect;connect(i,:)];
%     else
%         newbran{i}=[];
%     end
% end
% id = cellfun('length',newbran);newbran(id==0)=[];    
numbran = length(newbran);
middlebran = zeros(numbran,1);
for i = 1:numbran
    element = newbran{i};
    lensbran = size(element,1);
    if (mod(lensbran,2)==1)%odd
        middle = floor(lensbran/2)+1;
        middlebran(i) = element(middle);
    else %even
        middle = floor(lensbran/2);
        middlebran(i)=element(middle);
    end
end
% direction = zeros(numbran,3);
% [x1,y1,z1] = ind2sub(datas,newconnect(:,1));
% [x2,y2,z2] = ind2sub(datas,newconnect(:,2));
% aa = abs([x1-x2,y1-y2,z1-z2]);
% direction(:,1:3) = aa./ sqrt(sum(abs(aa).^2,2));
%%
newbran_thickness=cell(1);
newconnect_thickness = [];
%%
index =1;
clear thickness
for i = 1: length(middlebran)
    [x0,y0,z0] = ind2sub(datas,middlebran(i,1));
    width = 15;
%     if x0>width && x0<datas(1)-width &&y0>width && y0<datas(2)-width &&z0>width && z0<datas(3)-width
%         fprintf('x0 y0 z0 %d %d %d\n',x0,y0,z0)
        w = 1;
        temp = 0;
        while temp <1
            w = w+1;
            region = cont(x0-w:x0+w,y0-w:y0+w,z0-w:z0+w);
            temp = length(find(region == 1));
        end
%         fprintf('iteration %d w %d\n',i,w)
    newbran_thickness{index} = newbran{i};
    newconnect_thickness(index,:) = newconnect(i,:);
    thickness(index,1) = middlebran(i,1);
    thickness(index,2) = w;%pi*radius^2;
    index = index+1;
%     end
end       