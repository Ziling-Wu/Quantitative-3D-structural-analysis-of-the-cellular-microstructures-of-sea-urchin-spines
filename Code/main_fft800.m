clearvars
addpath('/Volumes/GoogleDrive/My Drive/Research/Dr. Li collaboration/Sea Urchin/20190301_Center/fft800/segmentation800')
%% load images
tic
startpage = 0;
endpage = 814;
style = '.tif';
for i = startpage:endpage
    if i<10
        im = imread(fullfile( ['recon_proj_0005_000' num2str(i) style]))*1e4;
    elseif i>=10 && i<100
        im = imread(fullfile( ['recon_proj_0005_00' num2str(i) style]))*1e4;
    elseif i>=100
        im = imread(fullfile( ['recon_proj_0005_0' num2str(i) style]))*1e4;
    end
    im(im==255)=1;
    seg(:,:,i+1) = im;
end
toc
%%
tic
skel_ML = Skeleton3D(seg);
toc
%% Trimming
tic
datas = size(skel_ML);
global datas
trimnode = node_identification(skel_ML);
trimskel = skel_ML;
% trim the branches and renode 
criteria = 2;
criteria_new = 1;
while criteria ~= criteria_new
    fprintf('1')
    criteria = size(trimnode,1);
    trimskel = branch_trim(trimnode, trimskel);
%     trimskel = Skeleton3D(trimskel);
    trimnode = node_identification(trimskel);
    criteria_new = size(trimnode,1);
end
toc
%%
save('trimnode_fft800.mat','trimnode')
save('skel_fft800.mat','skel_ML')
save('trimskel_fft800.mat','newtrimskel')
%%
tic
[bran,connect, node] = branch_sort2(trimskel,trimnode);
toc
%%
save('node.mat','node')
save('connect.mat','connect')
save('bran.mat','bran')
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
%% Fourier trnasform
fftnode = abs((ifftshift(fftn(fftshift(newskel)))));%*(0.65^3);
fftnode = fftnode - min(fftnode(:));
fftnode = fftnode./max(fftnode(:));
%%
save('fftnode.mat','fftnode')
%%
figure
thr = 0.02;
cor = find(fftnode>thr);
[a,b,value] = find(fftnode>thr);
[cor1,cor2,cor3] = ind2sub(datas,cor); % find nonzero component
scatter3(cor2(:,1),cor1(:,1),cor3(:,1),10,value,'filled');
colormap gray
axis image
set(gcf,'color','white')
%%

for i =1:815
    im = fftnode(:,:,i);
    imwrite(im,[num2str(i),'.tif']);
end