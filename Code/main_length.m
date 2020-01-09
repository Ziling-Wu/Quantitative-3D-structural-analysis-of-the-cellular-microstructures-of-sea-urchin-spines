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
innerbran = cell(1);
innerconnect = [];
for i=1:size(bran,2)
    bic(i)=isboundary(bran{i},30); %if 1, branch is at the boundary.
    if ~bic(i)   
        innerbran{index}=bran{i};
      	innerconnect=[innerconnect;connect(i,:)];
      	index=index+1;
	end
end  
clear bran
clear connect
connect = innerconnect;
bran = innerbran;

%% bran length: branlen, and branch distance: brandis
branlen = [];
brandis = [];
for i = 1: length(bran)
    tmpbran= bran{i};
    dis = 0;
    for j =1:length(tmpbran)-1
        dis = vecdist(tmpbran(j),tmpbran(j+1))+dis;
    end
    startDis = [vecdist(connect(i,1),tmpbran(1)),vecdist(connect(i,1),tmpbran(end))];
    EndDis = [vecdist(connect(i,2),tmpbran(1)),vecdist(connect(i,2),tmpbran(end))];
% [xs, ~] = sort(startDis);
    sDis = min(startDis);
    eDis = min(EndDis);
    branlen = [branlen;dis+sDis+eDis];
%     branlen(i) = length(tmpbran);
%     Branch distance
    brandis = [brandis;vecdist(connect(i,1),connect(i,2))];
end
branlen = branlen*0.65;
brandis = brandis *0.65;
ratio = branlen./brandis;
%% delete repeated branch
[bran,connect] = DeleteRepeatedBran(connect,bran,ratio);
%% Calcualte again: bran length: branlen, and branch distance: brandis
branlen = [];
brandis = [];1
for i = 1: length(bran)
    tmpbran= bran{i};
    dis = 0;
    for j =1:length(tmpbran)-1
        dis = vecdist(tmpbran(j),tmpbran(j+1))+dis;
    end
    startDis = [vecdist(connect(i,1),tmpbran(1)),vecdist(connect(i,1),tmpbran(end))];
    EndDis = [vecdist(connect(i,2),tmpbran(1)),vecdist(connect(i,2),tmpbran(end))];
% [xs, ~] = sort(startDis);
    sDis = min(startDis);
    eDis = min(EndDis);
    branlen = [branlen;dis+sDis+eDis];
%     branlen(i) = length(tmpbran);
%     Branch distance
    brandis = [brandis;vecdist(connect(i,1),connect(i,2))];
end
branlen = branlen*0.65;
brandis = brandis *0.65;
ratio = branlen./brandis;
%% Plot: the distance and length distribuiton
clf
hdc= histofrequency(brandis,0:5:100);
% the length distribution
hlc= histofrequency(branlen,0:5:100);
plot(hdc(:,1),hdc(:,2),hlc(:,1),hlc(:,2));
legend('l_{ij}','l_{o,ij}')
xlabel('Branch length (\mum)');ylabel('Frequency')
axis square
xticks([0 20 40 60 80 100])
xlim([0,60])
ylim([0,0.35])
set(gcf,'color','white')
set(gca,'fontsize', 16);
%  saveas(gcf,'branch distance-length comparsion','epsc')
% saveas(gcf,'branch distance-length comparsion','pdf')
% saveas(gcf,'branch distance-length comparsion','png')
%% ratio
ratio = branlen./brandis;
hratio= histofrequency(ratio,0.95:0.03:max(ratio));

plot(hratio(:,1),hratio(:,2),'LineWidth',0.75)
xlim([0.9,2])
ylim([0,0.23])
xlabel('Branch length/branch distance');ylabel('Frequency')
axis square
yticks([0:0.05:0.2])
set(gcf,'color','white')
% saveas(gcf,'branch length distance ratio','epsc')
% saveas(gcf,'branch length distance ratio','pdf')
% saveas(gcf,'branch length distance ratio','png')

%% plot the branches based on ratio
clf
color=jet(30);
%
clf
cor = [];
for i=1:30
    group=((1+(i-1)*0.0167)<ratio).*(ratio<(1+i*0.0167));
    group=find(group);
    temp=[];
    if ~isempty(group)
        for j=1:size(group,1)
            temp=[temp;bran{group(j)}];
            cor = [cor;group(j)];
        end

        if ~isempty(temp)
        PlotScatter(temp,5,color(i,:),0)
        hold on
        end
    end
end
% group=(ratio>(1+30*0.0167));
% group=find(group);
% for j=1:size(group,1)
%             temp=[temp;bran{group(j)}];
%             cor = [cor;group(j)];
% end
% PlotScatter(temp,5,color(31,:),0)
% hold on
colormap(color) ;
cbh = colorbar ; %Create Colorbar
% set(cbh,'XTickLabel',{'0','10','20','30','40','50','60'})
cbh.TickLabels = num2cell(linspace(1,1.5,11));
PlotScatter(connect(cor(:),1),10,'r',0)
hold on
PlotScatter(connect(cor(:),2),10,'r',0)

% colormap jet(100)

FigureFormat
% saveas(gcf,'ratio length_center','epsc');
% saveas(gcf,'ratio length_center','pdf');
% saveas(gcf,'ratio length_center','png');
%% plot the branches based on different branch: this will take some time
clf
for i = 1:size(bran,2) 
    temp = bran{i};
    PlotScatter(temp,5,[1-i/size(bran,2),i/size(bran,2),1-i/size(bran,2)],0)
    hold on 
end

PlotScatter(connect(:,1),8,'r',0)
hold on
PlotScatter(connect(:,2),8,'r',0)
axis off
colormap jet(100)
FigureFormat

%% plot the branches based on the length center part
clf
color1=jet(30);
%
clf
cor = [];
for i=1:30
    group=((i-1)*2<branlen).*(branlen<i*2);
    group=find(group);
    temp=[];
    if ~isempty(group)
        for j=1:size(group,1)
            temp=[temp;bran{group(j)}];
            cor = [cor;group(j)];
        end

        if ~isempty(temp)
        PlotScatter(temp,5,color1(i,:),0)
        hold on
        end
    end
end

PlotScatter(connect(cor(:),1),10,'r',0)
hold on
PlotScatter(connect(cor(:),2),10,'r',0)
%  cmap = colormap(color) ; %Create Colormap
%  cbh = colorbar ; %Create Colorbar
%  cbh.Ticks = linspace(0, 60,10) ; %Create 8 ticks from zero to 1
%  cbh.TickLabels = num2cell(0:10:60);
% colormap jet(100)
cmap = colormap(color1) ;
cbh = colorbar ; %Create Colorbar
% set(cbh,'XTickLabel',{'0','10','20','30','40','50','60'})
cbh.TickLabels = num2cell(0:6:60);
     
FigureFormat
% colorbar
% caxis([10,30*2.5])
% saveas(gcf,'branch length_center','epsc');
% saveas(gcf,'branch length_center','pdf');
% saveas(gcf,'branch length_center','png');
%% center
figure
cor = [];
for i=1:7
% for i = 20:30
    group=((i-1)*2<branlen).*(branlen<i*2);
    group=find(group);
    temp=[];
    if ~isempty(group)
        for j=1:length(group)
            temp=[temp;bran{group(j)}];
            cor = [cor;group(j)];
        end  
        if ~isempty(temp)
        PlotScatter(temp,5,color1(i,:),0)
        hold on
        end
    end
end

PlotScatter(connect(cor(:),1),10,'r',0)
hold on
PlotScatter(connect(cor(:),2),10,'r',0)
% caxis([10,30*2.5])
FigureFormat
axis([1 308 1 386 1 386])
saveas(gcf,'branch length_center_short_branches','epsc');
saveas(gcf,'branch length_center_short_branches','pdf');
saveas(gcf,'branch length_center_short_branches','png');
%% Correlation: branch orientation vs branch length analysis_center

clf
record = [];
for i=1:30
    group=((i-1)*2<branlen).*(branlen<i*2);
    group=find(group);
    temp=[];
    if ~isempty(group)
        for j=1:length(group)    
            [theta,phi]=myangle(connect(group(j),1),connect(group(j),2));
            temp(j,:)=real([theta,phi]);
        end
        if ~isempty(temp)
            record = [record;temp];
            p=polarplot(deg2rad(temp(:,1)),abs((temp(:,2))),'o','MarkerSize',6);
            p.Color = color1(i,:);
            hold on
        end
    end    
 
end
% thetalim([-90 90]);
 %colorbar
% caxis([0,22*2.5])
a=gca;
 a.LineWidth=2;
 a.FontSize=12;
 set(gcf,'color','white')
% saveas(gcf,'Center region branch orientation_without colorbar','epsc')
%  saveas(gcf,'Center region branch orientation_without colorbar','pdf')
%  saveas(gcf,'Center region branch orientation_without colorbar','png')
 
 
% %%
% clear temp
% index=1;
% for i =1:length(bran)
%      [theta,phi]=myangle(connect(i,1),connect(i,2));
%       temp(i,:)=real([theta,phi]);
%       if temp(i,1)>90
%           tmpbran = [tmpbran;bran{i}];
%           index = index+1;
%       end
% end
% PlotScatter(tmpbran,15,'g')
% hold on
% PlotScatter(trimskel,10,'r')
%%
% histogram(temp(:,1))