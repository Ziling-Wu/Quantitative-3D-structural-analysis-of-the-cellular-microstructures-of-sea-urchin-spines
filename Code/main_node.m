%% Branch length
clearvars
load('c5.mat');
load('b5.mat');
load('n5.mat');
load segML.mat
img = segML;
node = n5;
connect = c5;
bran = b5;
load trimskel
datas = size(trimskel);
global datas
%% is near boundary analysis 
index = 1;
clear innerNode
for i=1:size(node,1)
    nic(i)=isboundary(node(i,1),15); %if 1, branch is at the boundary.
    if ~nic(i)   
        innerNode(index,:) = node(i,:);
      	index=index+1;
	end
end  
clear node
node =innerNode;
%% Look at node type at first. 
hlc= histofrequency(node(:,2),3:9);
x=hlc(:,1)-0.5;
bar(x,hlc(:,2))
colormap jet
set(gcf,'color','white')
%% interbranch angle analysis 
nc = node;
cc = connect;
bc =bran;
clear angle_c
for i=1:size(nc,1)
    temp=cc(cc(:,1)==nc(i,1),2);
    temp=[temp;cc(cc(:,2)==nc(i,1),1)];
    for j=1:length(temp)
        angle_c(i,j)=interangle(nc(i,1),temp(j),temp(mod(j,length(temp))+1))*180/pi;
    end
end
%% 3 node angle 
% all angle 
clf
clear temp hlc

edge=[20:10:180]-2.5;
angle3N=abs(angle_c((angle_c(:,4)==0)&(angle_c(:,3)~=0),1:3));
temp=mean(angle3N,2);
hlcm= histcounts(real(temp),edge);
temp=min(angle3N,[],2); % smallest angle;
hlcs= histcounts(real(temp),edge);
temp=sort(angle3N,2);
temp=temp(:,2); % middle angle 
hlcmi= histcounts(real(temp),edge);
temp=max(angle3N,[],2); % largest angle; 
hlcl= histcounts(real(temp),edge);
x=edge(1:end-1)+2.5;
plot(x,hlcm,x,hlcs,x,hlcmi,x,hlcl,'LineWidth',0.75);
legend('\gamma','\gamma_1','\gamma_2','\gamma_3')
xlabel('Interbranch angles (degrees)');ylabel('Counts')
xticks([20,60,100,140,180])
yticks([0:100:450])
ylim([0 50])
set(gca,'FontSize',7);
set(0,'DefaultTextFontName','Ariel')
set(gcf,'color','white')
% saveas(gcf,'3-node interbranch angle','epsc')
% saveas(gcf,'3-node interbranch angle','png')
% saveas(gcf,'3-node interbranch angle','pdf')
%%
temp=mean(angle3N,2);
mean(angle3N(:))
std(temp)
%% 4 node angle 
clf
clear temp hlc hle
edge=[20:10:180]-5;
angle4N=abs(angle_c((angle_c(:,4)~=0)&(angle_c(:,5)==0),1:4));
temp=mean(angle4N,2);
hlcm= histcounts(temp,edge);
%
temp=max(angle4N,[],2); % largest angle; 
hlcl= histcounts(real(temp),edge);
temp=min(angle4N,[],2); % largest angle; 
hlcs= histcounts(real(temp),edge);
x=edge(1:end-1)+5;

plot(x,hlcm,x,hlcs,x,hlcl,'LineWidth',0.75);
legend('\gamma','\gamma_{min}','\gamma_{max}')
xlabel('Interbranch angles (degrees)');ylabel('Counts')
set(gcf,'color','white')
set(gca,'FontSize',7);
set(0,'DefaultTextFontName','Ariel')
xticks([20,60,100,140,180])
yticks([0:200:1400])
ylim([0 1400])
% saveas(gcf,'4-node interbranch angle','epsc')
% saveas(gcf,'4-node interbranch angle','png')
% saveas(gcf,'4-node interbranch angle','pdf')
%%
temp=mean(angle4N,2);

mean(angle4N(:))
std(temp)
%% node type
clf

colorMap=[1,0.6,0.1;1,1,0;1,0,1;0,0,1];
for i=1:4
PlotScatter(nc(nc(:,2)==i+2,1),25,colorMap(i,:),0)
hold on
end
hold on
temp=[];
for i=1:size(bc,2)
    temp=[temp;bc{i}];
end
hold on
PlotScatter(temp,1.5,[0.5,0.5,0.5],0)
FigureFormat
colormap spring(4)
%colorbar 
% saveas(gcf,'nodetype_center_thinner','epsc');
% saveas(gcf,'nodetype_center_thinner','pdf');
% saveas(gcf,'nodetype_center_thinner','png');
%% plot based on different branch
clf
temp =[];
for i = 1:size(bc,2) 
    hold on
    temp = [temp;bc{i}];
    
end
PlotScatter(temp,1.5,[0.502,0.502,0.502],0)
    hold on 
PlotScatter(n5(:,1),25,'r',0)
FigureFormat
axis image
% saveas(gcf,'network','epsc');
% saveas(gcf,'network','png');
% saveas(gcf,'network','pdf');
%% node type
clf
hlc= histofrequency(nc,3:1:7);
x=hlc(:,1)-0.5;
bar(x,hlc(:,2))
colormap jet
set(gcf,'color','white')
xlabel('Node types');ylabel('Frequency')
ylim([0,0.6])
% saveas(gcf,'nodetype','epsc')
% saveas(gcf,'nodetype','png')
% saveas(gcf,'nodetype','pdf')
%% the equal distance orientation
tic
clear oric
clear plane
n3c=find(nc(:,2)==3);
clear temp
for i=1:size(n3c,1)
    temp=cc(cc(:,1)==nc(n3c(i),1),2);
    temp=[temp;cc(cc(:,2)==nc(n3c(i)),1)];
    [oric(i,:),plane(i,:)] =ori3(nc(n3c(i),1),temp(1),temp(2),temp(3));
end
toc
%% center
% load oric
figure
 clear temp
 temp(:,2)=atand(sqrt(oric(:,3).^2+oric(:,2).^2)./abs(oric(:,1)));% phi
% temp(:,2)=acosd(sqrt(oric(:,3).^2+oric(:,2).^2));% phi
 for i =1:size(oric,1)
    % temp(:,1)=atan(oric(:,2)./oric(:,3)); 
    theta= abs(atan(abs(oric(i,2))/abs(oric(i,3))))*180/pi;
    if oric(i,1)>0
        temp(i,2) = temp(i,2);
        y = oric(i,2);
        z = oric(i,3);
    else
        temp(i,2) = 180-temp(i,2);
        y = oric(i,2);
        z = oric(i,3);
    end
    if (y>0)&&(z>0)
        theta = 360-theta;
    elseif (y>0)&&(z<0)
        theta = 180+theta;
    elseif (y<0)&&(z<0)
        theta = 180-theta;
    else
        theta = theta;
    end
    temp(i,1) = theta;
 end
clf
 p=polarplot(deg2rad(temp(:,1)),temp(:,2),'o','MarkerSize',6,'LineWidth',0.5);
rticks([0:30:180])
 set(gcf,'color','white')
 saveas(gcf,'Center node orientation','epsc')
saveas(gcf,'Center node orientation','pdf')
saveas(gcf,'Center node orientation','png')
%%

figure
nbins = -0.009:0.05:1;
hlc= histofrequency(abs(cos(deg2rad(plane))),nbins);
plot(hlc(:,1),hlc(:,2),'LineWidth',1)
x=hlc(:,1);
% bar(x,hlc(:,2))
colormap jet
set(gcf,'color','white')
xlabel('Planarity angles');ylabel('Frequency')
xticks([0:0.2:1])
% yticks([0,0.05,0.1,0.15,0.2])
xlim([0,1])
% histogram(plane,nbins)
ylim([0,0.5])
% yticks([0:0.1:0.4])
set(gcf,'color','white')
 saveas(gcf,'Planarity','epsc')
saveas(gcf,'Planarity','pdf')
saveas(gcf,'Planarity','png')