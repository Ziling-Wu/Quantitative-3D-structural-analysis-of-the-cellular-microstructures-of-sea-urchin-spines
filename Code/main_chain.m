%% Branch length
clearvars
load('c5.mat');
load('b5.mat');
load('n5.mat');
node = n5;
connect = c5;
bran = b5;
load trimskel
datas = size(trimskel);
global datas
%% is near boundary analysis 
index = 1;
innerconnect = [];
innerbran = cell(1);
for i=1:size(connect,1)
    bic(i)=isboundary(connect(i,:),10); %if 1, branch is at the boundary.
    if ~bic(i)   
      	innerconnect=[innerconnect;connect(i,:)];
        innerbran{index} = bran{i};
      	index=index+1;
	end
end  
clear connect
clear bran
bran = innerbran;
connect = innerconnect;
%%
branOrien = zeros(length(bran),1);
for i=1:length(bran)
    tempbran=bran{i};   
    [theta,phi]=myangle(tempbran(1),tempbran(end));   
    branOrien(i) = phi;
end
%%
node = unique(connect);
% startnode
[x,y,z] =  ind2sub(datas,node(:,1));
startnode = [];
for i = 1:length(x)
    if (x(i)<40)&&(x(i)>10)
        startnode = [startnode;node(i,1)];
    end
end
%
allchain = cell(length(startnode),1);
for i = 1:size(startnode,1)
    chainnode = [];
    chainnode = chain(startnode(i),connect,chainnode,branOrien);
    allchain{i}=chainnode;
end
%  cut
for j = 1:size(allchain,1)
    chainnode = allchain{j,1};
    loc = find(chainnode==0);
    loc = [0;loc];
    longline = cell(length(loc)-1,1);
    for i = 1: length(loc)-1
        longline{i} = chainnode(loc(i)+1:loc(i+1)-1);
    end
    allchain{j,2} = longline;
end
% stitich
for j = 1:size(allchain,1)
    longline = allchain{j,2};
    for i = 2:length(longline)
        firstline = longline{i-1};
        temp = longline{i};
        cut = find(firstline == temp(1));
        stitch = [firstline(1:cut);temp(2:end)];
        longline{i} = stitch;
    end
    allchain{j,2} = longline;
end
% longest chain
for j = 1:size(allchain,1)
    longline = allchain{j,2};
    len = [];
    for i = 1:size(longline,1)
        len(i) = length(longline{i});
    end
    loc_line = find(len == max(len));
    len_line = max(len);
    allchain{j,3} = loc_line;
    allchain{j,4} = len_line;
end
%%
len_line = [];
for i = 1:size(allchain,1)
    len_line = [len_line,allchain{i,4}];
end
loc_line = find(len_line >=10);
%%
clf
% tmpbran = [];
% for i =1:length(bran)
%     tmpbran = [tmpbran;bran{i}];
% end
PlotScatter(trimskel,5,[1,1,1]-0.3,0);
hold on 
tmpbran = [];
for k = 1:length(loc_line)
    branNumber = [];
    loc_chain = allchain{loc_line(k),3};
    temp = allchain{loc_line(k),2};
    for i = 1:length(loc_chain)
        chaindots = temp{loc_chain(i)};
        for j = 1:size(chaindots,1)-1
            c = findbran(chaindots(j),chaindots(j+1),connect);
%             clear c1 c2
%             c1 = find((connect(:,1)==chaindots(j))&(connect(:,2)==chaindots(j+1)));
%             c2 = find((connect(:,2)==chaindots(j))&(connect(:,1)==chaindots(j+1)));
%             c = [c1;c2];
            t = bran{c};
            tmpbran = [tmpbran;t];
        end
        index = index+1;
        PlotScatter(temp{loc_chain(i)},50,'r',0)
        hold on
    end
end
PlotScatter(tmpbran,10,'b',0)
axis on
FigureFormat(90,90)
saveas(gcf,'chain','epsc');
saveas(gcf,'chain','pdf');
saveas(gcf,'chain','png');
% saveas(gcf,'long_chain','epsc');
% saveas(gcf,'long_chain','pdf');
% saveas(gcf,'long_chain','png');
 %%
clear offangle
clear chaindirec
chainleng = zeros(length(loc_line),1);
ind = 1;
chaindirec = cell(1);
 for k = 1:length(loc_line)
    loc_chain = allchain{loc_line(k),3};
    temp = allchain{loc_line(k),2};
    len = allchain{loc_line(k),4};
    for i = 1:length(loc_chain)
        chaindots = temp{loc_chain(i)};
        for j = 1:size(chaindots,1)-1
            c = findbran(chaindots(j),chaindots(j+1),connect);
            t = bran{c};
            dis = 0;
            for kk =1:length(t)-1
                dis = vecdist(t(kk),t(kk+1))+dis;
            end
            startDis = [vecdist(chaindots(j),t(1)),vecdist(chaindots(j),t(end))];
            EndDis = [vecdist(chaindots(j+1),t(1)),vecdist(chaindots(j+1),t(end))];
% [xs, ~] = sort(startDis);
            sDis = min(startDis);
            eDis = min(EndDis);
            branlen = dis+sDis+eDis;
            [cx1,cy1,cz1] = ind2sub(datas,t(1));
            [cx2,cy2,cz2] = ind2sub(datas,t(end));
            chainleng(ind,j) = branlen;
            chaindirec{ind,j} = [cx2-cx1,cy2-cy1,cz2-cz1];
            offangle(ind,j) = branOrien(c);
        end
    end
    ind = ind+1;
end
%%
figure
stas = nonzeros(chainleng*0.65);
mean(stas)
histogram(stas(stas<=50));
xlim([0,45])
ylim([0,160])
xticks([0:10:50])
yticks([0:20:160])
set(gcf,'color','white')
saveas(gcf,'chainleng','epsc');
saveas(gcf,'chainleng','pdf');
saveas(gcf,'chainleng','png');
%% direction
chainAngle = zeros(length(loc_line),1);
ind = 1;
for i = 1:size(chaindirec,1)
    for j = 1:size(chaindirec,2)-1
        x1 = cell2mat(chaindirec(i,j));
        x2 = cell2mat(chaindirec(i,j+1));
        if ~isempty(x2)
            temp = acosd(dot(x1,x2)/(norm(x1)*norm(x2)));
            if temp >= 90
                chainAngle(i,j) = temp;
            else
                chainAngle(i,j) = 180-temp;
            end
        else
            continue
        end
    end
end
%%
figure
nbins = 90:5:180;
stas = nonzeros(chainAngle);
histogram(stas,nbins)
ylim([0,300])
set(gcf,'color','white')
xticks([80:10:180])
yticks([0:50:300])
saveas(gcf,'InterchainAngle','epsc');
saveas(gcf,'InterchainAngle','pdf');
saveas(gcf,'InterchainAngle','png');
%% off-L direction
figure
nbins = 0:5:50;
histogram(offangle(:),nbins)
ylim([0,300])
set(gcf,'color','white')
xticks([0:10:50])
yticks([0:50:300])
saveas(gcf,'offangle','epsc');
saveas(gcf,'offangle','pdf');
saveas(gcf,'offangle','png');
%% branch thickness
tic
load cont
chainThick = cell(length(loc_line),1);
index = 1;
for k = 1:length(loc_line)
    loc_chain = allchain{loc_line(k),3};
    tempCell = allchain{loc_line(k),2};
    for i = 1:length(loc_chain)
        chaindots = tempCell{loc_chain(i)};
        for j = 1:size(chaindots,1)-1
            c = findbran(chaindots(j),chaindots(j+1),connect);
            tmpbran = bran{c};
            % find general thickness
            [x0,y0,z0] = ind2sub(datas,tmpbran(1));
            w = 1;temp = 0;
            while temp <5
                w = w+1;
                region = cont(x0-w:x0+w,y0-w:y0+w,z0-w:z0+w);
                temp = length(find(region == 1));
            end  
            % find thickness
            [x,y,z]=meshgrid(1:datas(1),1:datas(2),1:datas(3));
            [x1,y1,z1] = ind2sub(datas,tmpbran(1));
            [x2,y2,z2] = ind2sub(datas,tmpbran(end));
            r1 = [x2-x1,y2-y1,z2-z1];
            l = length(tmpbran);
            if l<5
                continue
            else
                ind = 1;
                thick = [];
                for jj = 2:1:length(tmpbran)-1
                    [x0,y0,z0] = ind2sub(datas,tmpbran(jj));
                    tempsphere = cont.*(((x-y0).^2+(y-x0).^2+(z-z0).^2)<=20^2);
                    cor = find(tempsphere);
                    if isempty(cor)
                        continue
                    else
                        [xr,yr,zr] = ind2sub(datas,cor);
                        c = [xr-x0,yr-y0,zr-z0];
                        clear angle1
                        angle1 = abs(acos(c*r1'/norm(r1)./vecnorm(c,2,2))*180/pi)-90;
                        loc1 = find(abs(angle1)<2);% == min(abs(angle)));
                        loc = loc1;
                        if isempty(loc)
                            fprintf('here!\n');
                            continue;
                        end
                        dis = sqrt((x0-xr(loc)).^2+(y0-yr(loc)).^2+(z0-zr(loc)).^2);
                        if sum(dis<w+11) == 0
                            thick(ind) = mean(dis);
                        else
                            thick(ind) = mean(dis((dis<w+11)));
                        end
                        ind = ind+1;
                    end
                end
                chainThick{index,j} = thick;
            end
        end
    end
    index = index+1;
end
toc
%%
save('chainThick.mat','chainThick')
%%
chainMinthickness = zeros(size(chainThick));
for i = 1:size(chainThick,1)
    for j = 1:size(chainThick,2)
        tempRow = chainThick{i,j};
        if ~isempty(tempRow)
            chainMinthickness(i,j) = min(tempRow)*0.65;
        end
    end
end
figure
allChainthickness = chainMinthickness(chainMinthickness~=0);
mean(allChainthickness)
std(allChainthickness)
nbins = [1.5:1:10.5];
histogram(allChainthickness,nbins);
xlim([2,11])
ylim([0,350])
xticks([3:2:10])
% yticks([0:20:150])
% saveas(gcf,'chainthick','epsc');
% saveas(gcf,'chainthick','pdf');
% saveas(gcf,'chainthick','png');
%%

allChainlength = chainleng(chainMinthickness~=0)*0.65;
mean(allChainlength)
std(allChainlength)
scatter(allChainlength,allChainthickness)
xlim([0,40])
% ylim([0,280])
xticks([0:10:50])
yticks([0:2:10])
xlabel('Branch length')
ylabel('Branch thickness')
% saveas(gcf,'chainThickLength','epsc');
% saveas(gcf,'chainThickLength','pdf');
% saveas(gcf,'chainThickLength','png');
%%
save('allChainlength.mat','allChainlength');
save('allChainthickness.mat','allChainthickness');
save('chainThick.mat','chainThick');