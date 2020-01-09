
%% load images
tic
startpage = 0;
endpage = 307;
style = '.tif';
for i = startpage:endpage
    if i<10
        im = imread(fullfile( ['000' num2str(i) style]))*1e4;
    elseif i>=10 && i<100
        im = imread(fullfile( ['00' num2str(i) style]))*1e4;
    elseif i>=100
        im = imread(fullfile( ['0' num2str(i) style]))*1e4;
    end
    im(im==255)=1;
    segML(:,:,i+1) = im;
end
toc
%%
save('segML.mat','segML')
%%
% %% skeltonize
skel_ML = Skeleton3D(segML);
%%
save('skel_ML.mat','skel_ML')
% PlotScatter(trimskel,2,'r')

%% Trimming
datas = size(skel_ML);
global datas
trimnode = node_identification(skel_ML);
trimskel = skel_ML;
%% trim the branches and renode 
tic
criteria = 2;
criteria_new = 1;
while criteria ~= criteria_new
    fprintf('1\n')
    criteria = size(trimnode,1);
    trimskel = branch_trim(trimnode, trimskel);
%     trimskel = Skeleton3D(trimskel);
    trimnode = node_identification(trimskel);
    criteria_new = size(trimnode,1);
end
toc

%% identify the nodes and register branches
datas = size(trimskel);
global datas
[bran,connect, node] = branch_sort2(trimskel,trimnode);
bran = branch_order(bran,connect,datas);
save('bran.mat','bran');
save('connect.mat','connect');
save('node.mat','node');
%% merge connected nodes
clearvars
load trimskel
datas = size(trimskel);
global datas
load node
load connect
load bran
[c2,n2,b2] = merge_connectednodes(connect,node,bran);

%% merge nearby nodes
[c3,n3,b3] = merge_Allnear(c2,n2,b2,10);
%% delete branches
[c4,n4,b4] = merge_double(c3,n3,b3);
%% merge nearby nodes
[c5,n5,b5] = merge_Allnear(c4,n4,b4,10);
%%
nodef = [];
for i = 1:size(n5,1)
    if n5(i,2)>2
        nodef = [nodef;n5(i,:)];
    end
end
%
n5 = nodef;
%%
save('trimskel.mat','trimskel');

%%
save('b5.mat','b5');
save('c5.mat','c5');
save('n5.mat','n5');
%%
save('b2.mat','b2');
save('c2.mat','c2');
save('n2.mat','n2');
% %% bran length;
% clearvars
% load('c5.mat');
% load('b5.mat');
% load('n5.mat');
% datas = [386,386,308];
% global datas
% branlen = [];
% brandis = [];
% for i = 1: length(b5)
%     
% % Branch length
%     tmpbran= b5{i};
%     dis = 0;
%     [x,y,z] = ind2sub(datas,tmpbran);
%     if ((x<20)|(x>360)|(y<20)|(y>360)|(z<20)|(z>280))
%         continue
%     else
%     for j =1:length(tmpbran)-1
%         dis = vecdist(tmpbran(j),tmpbran(j+1))+dis;
%     end
%     startDis = [vecdist(c5(i,1),tmpbran(1)),vecdist(c5(i,1),tmpbran(end))];
%     EndDis = [vecdist(c5(i,2),tmpbran(1)),vecdist(c5(i,2),tmpbran(end))];
% % [xs, ~] = sort(startDis);
%     sDis = min(startDis);
%     eDis = min(EndDis);
%     branlen = [branlen;dis+sDis+eDis];
% %     branlen(i) = length(tmpbran);
% %     Branch distance
%     [y1,x1,z1] = ind2sub(datas,c5(i,1));
%     [y2,x2,z2] = ind2sub(datas,c5(i,2));
% %     connectivity = [x1,y1,z1,x2,y2,z2];
%     brandis = [brandis;vecdist(c5(i,1),c5(i,2))];
%     end
% end
% % histogram(branlen)
% % hold on
% histogram(brandis,'FaceColor','r')
% hold on
% histogram(branlen,'FaceColor','g')
% %%
% ratio = branlen./brandis;
% %%
% [bran,connect] = DeleteRepeatedBran(c5,b5,ratio);


% %% 
% outSkel = uint8(zeros(datas));
% tmpbran =[];
% for i =1:length(bran)
%      tmpbran = [bran{i};tmpbran];
% end
% outSkel(tmpbran)=1;
% outSkel(node(:,1))=10;
% %%
% PlotScatter(trimnode(:,1),2,'r')
% %%
% 
% for i =1:308
%     im = outSkel(:,:,i);
%     imwrite(im,[num2str(i),'.tif']);
% end
