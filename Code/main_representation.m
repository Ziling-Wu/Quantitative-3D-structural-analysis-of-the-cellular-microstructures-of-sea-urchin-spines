%% merge connected nodes
clearvars
load node
load connect
load bran
load trimskel_fft800
datas = size(newtrimskel);
global datas
[c2,n2,b2] = merge_connectednodes(connect,node,bran);
%%
newskel = zeros(datas);
newskel(n2(:,1))=1;
% fftnode = ((ifftshift(fftn(fftshift(newskel)))));%*(0.65^3);
% %%
% save('fftnode.mat','fftnode')
%%
load fftnode
%%
fftnode2 = imresize3(fftnode,0.5);
%%
subplot(121);imagesc(abs(fftnode2(:,:,1)));axis image
subplot(122);imagesc(abs(fftnode(:,:,1)));axis image
%% Control node periodity by delete these small value frequency
figure
mag=abs(fftnode2);
histogram(mag(:))
%% 
fftnode3 = fftnode2;

fftnode3(abs(fftnode2(:))<=prctile(mag(:),99.999))=0;
inode = ifftshift(ifftn(fftshift((fftnode3))));%
%%
inode = inode - min(inode(:));
inode = inode./max(inode(:));
inodeP = find(inode>0.65);
datas = size(inode);
global datas
iskel = zeros(datas);
iskel(inodeP)=1;
figure
% subplot(122)
PlotScatter(iskel)
FigureFormat

% axis([1 200 1 200 1 200])
title('Reconstructed node Distribution,99.995%')
%%
subplot(121)
PlotScatter(newskel)
FigureFormat
title('Original node Distribution')
% axis([1 200 1 200 1 200])
%% Establish connectivity


inode = ifftshift(ifftn(fftshift((fftnode2))));%
inode = inode - min(inode(:));
inode = inode./max(inode(:));
inodeP = find(inode>0.6);
datas = size(inode);
global datas
iskel = zeros(datas);
iskel(inodeP)=1;
figure
% subplot(122)
PlotScatter(iskel)
FigureFormat
%%
% load n5
% load trimskel
% datas = size(trimskel);
% global datas
% newskel = zeros(size(trimskel));
% newskel(n5(:,1))=1;
cropSize = 99;star=100;
cropSkel = newskel(star:star+cropSize,star:star+cropSize,star:star+cropSize);
cropSkel_node = find(cropSkel);
% D = graph(cropSkel_node);
%% establish all connections
clear iconnect iconnect_cor
tmpConnect=[];
n_node = length(cropSkel_node);
for i = 1:n_node-1
    tmpConnect = [tmpConnect,repelem(i,n_node-i)];
end
iconnect(:,1)= tmpConnect';
for i = 1:n_node-1
    iconnect(iconnect(:,1)==i,2) = i+1:n_node;
end

%
iconnect_cor(:,1) = cropSkel_node(iconnect(:,1));  
iconnect_cor(:,2) = cropSkel_node(iconnect(:,2));  
%% Criterion
for j = 1:length(iconnect_cor)
    [x1,y1,z1] = ind2sub([100,100,100],iconnect_cor(j,1));
    [x2,y2,z2] = ind2sub([100,100,100],iconnect_cor(j,2));
    iconnect_cor(j,3)=sqrt((x1-x2).^2 + (y1-y2).^2 + (z1-z2).^2) ;
    iconnect_cor(j,4)=atand(sqrt((y1-y2).^2 + (z1-z2).^2)/abs(x1-x2));
end
idx=randperm(n_node);
n3 = sort(idx(1:round(n_node*0.6)));
n4 = sort(idx(round(n_node*0.6)+1:round(n_node*0.95)));
n5 = sort(idx(round(n_node*0.95)+1:end));
cropSkel_node(n3,2)=3;
cropSkel_node(n4,2)=4;
cropSkel_node(n5,2)=5;
%% 
bran_already = 0;
new_iconnect =[];
for i = 1:n_node-2 
    clear R
    loc = find(iconnect(:,1)==i);
    maxlen = max(iconnect_cor(loc(:),3));
    maxAngle = max(iconnect_cor(loc(:),4));
    iconnect_cor(loc(:),5) = 1-iconnect_cor(loc(:),3)/maxlen+1e-6;
    iconnect_cor(loc(:),6) = 1-iconnect_cor(loc(:),4)/maxAngle+1e-6;
    iconnect_cor(loc(:),7) = iconnect_cor(loc(:),5)+ iconnect_cor(loc(:),6)/5;
    iconnect_cor(loc(:),8) = iconnect_cor(loc(:),7)/sum(iconnect_cor(loc(:),7));
    if i >1
        bran_already = length(find(new_iconnect(:,2) == i));
    end
    if (bran_already>=cropSkel_node(i,2))||length(loc)<cropSkel_node(i,2)-bran_already
        continue
    elseif length(loc) == cropSkel_node(i,2)-bran_already
        R = loc;
    else
        R = sort(randsample(loc, cropSkel_node(i,2)-bran_already, true, iconnect_cor(loc(:),8)));
        while length(unique(R)) ~= (cropSkel_node(i,2)-bran_already)
            R = sort(randsample(loc, cropSkel_node(i,2)-bran_already, true, iconnect_cor(loc(:),8)))
        end
    end
    for j = 1:length(R)-1
        for k= 2:length(R)
            de_loc = find((iconnect(:,1)==R(j))&(iconnect(:,2)==R(k)))
            iconnect(de_loc,:)=[];
            iconnect_cor(de_loc,:)=[];
        end
    end
    new_iconnect = [new_iconnect;iconnect(R,1:2)];
end
if length(find(new_iconnect(:,2) == n_node))+length(find(new_iconnect(:,1) == n_node))<cropSkel_node(n_node,2)
    new_iconnect = [new_iconnect;iconnect(end,1:2)];
end
%
new_iconnect_cor = cropSkel_node(new_iconnect);
%%
clf
datas = size(cropSkel);
global datas
[x,y,z] = ind2sub(datas,new_iconnect_cor(:,1));
xyz1 = [x,y,z];
[x,y,z] = ind2sub(datas,new_iconnect_cor(:,2));
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
PlotScatter(cropSkel_node(:,1))
FigureFormat
%%
%%
load c5
clf
datas = size(cropSkel);
global datas
[x,y,z] = ind2sub(size(trimskel),c5(:,1));
xyz1 = [x,y,z];
[x,y,z] = ind2sub(size(trimskel),c5(:,2));
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
% PlotScatter(cropSkel_node(:,1))

FigureFormat
axis([star star+cropSize star star+cropSize star star+cropSize])

%%

%%
% figure
% [bran,connect, node] = branch_sort2(cropSkel,trimnode);
% [x,y,z] = ind2sub([50,50,50],connect(:,1));
% xyz1 = [x,y,z];
% [x,y,z] = ind2sub([50,50,50],iconnect(:,2));
% xyz2 = [x,y,z];
% for i = 1:size(xyz1,1)
%     xx = [xyz1(i,1),xyz2(i,1)];
%     yy= [xyz1(i,2),xyz2(i,2)];
%     zz= [xyz1(i,3),xyz2(i,3)];
%     plot3(xx,zz,yy)
% %     hold on
% %     plot3(,'r')
%     hold on
% end
% hold on 
% scatter3(xyz1(:,1),xyz1(:,3),xyz1(:,2))