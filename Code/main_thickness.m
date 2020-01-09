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
% %% Branch centeral thickness
% [thickness,newbran_thickness,newconnect_thickness] = branch_thickness(bran,connect,node,img);
% %%
% record = find(thickness(:,2)>=8);
% tmpbran =[];
% for i =1:length(record)
%      tmpbran = [tmpbran;bran{record(i)}];
% end
% PlotScatter(tmpbran,30,'b')
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
ratio = branlen./brandis;
%% delete repeated branch
[bran,connect] = DeleteRepeatedBran(connect,bran,ratio);

%% 
c3 = connect;
b3 = bran;
% calculate branch profile
load segML
% cont = edge3(segML,'approxcanny',0.6);
load cont
%%
index =1;
clear thickness
for i = 1: length(bran)
    tmpbran =bran{i};
    [x0,y0,z0] = ind2sub(datas,tmpbran(1));
    width = 15;
%     if x0>width && x0<datas(1)-width &&y0>width && y0<datas(2)-width &&z0>width && z0<datas(3)-width
%         fprintf('x0 y0 z0 %d %d %d\n',x0,y0,z0)
        w = 1;
        temp = 0;
        while temp <5
            w = w+1;
            region = cont(x0-w:x0+w,y0-w:y0+w,z0-w:z0+w);
            temp = length(find(region == 1));
        end
%         fprintf('iteration %d w %d\n',i,w)
    newbran_thickness{index} = i;
    newconnect_thickness(index,:) = connect(i,:);
    thickness(index,1) = tmpbran(1);
    thickness(index,2) = w;%pi*radius^2;
    index = index+1;
%     end
end
%%
tic
[x,y,z]=meshgrid(1:datas(1),1:datas(2),1:datas(3));
allthickness = cell(1);
index = 1;
for i = 1:size(c3,1)
    tmpbran = b3{i};
    [x1,y1,z1] = ind2sub(datas,tmpbran(1));
    [x2,y2,z2] = ind2sub(datas,tmpbran(end));
    r1 = [x2-x1,y2-y1,z2-z1];
    width = 10;
        l = length(tmpbran);
        if l<10
            continue
        else
%             if (ratio(i)>=1)&&(ratio(i)<=1.4)%not too curve branch
                [xx,yy,zz]= ind2sub(datas,tmpbran);
                ind = 1;
                thick=[];
                for j = 2:1:length(tmpbran)-1
                    [x0,y0,z0] = ind2sub(datas,tmpbran(j));
                    tempsphere = cont.*(((x-y0).^2+(y-x0).^2+(z-z0).^2)<=20^2);
                    cor = find(tempsphere);
                    if isempty(cor)
                        continue
                    else
                        [xr,yr,zr] = ind2sub(datas,cor);
                        c = [xr-x0,yr-y0,zr-z0];
                        clear angle1
                        angle1 = abs(acos(c*r1'/norm(r1)./vecnorm(c,2,2))*180/pi)-90;
%                         angle2 = abs(acos(c*r2'/norm(r2)./vecnorm(c,2,2))*180/pi)-90;
                        loc1 = find(abs(angle1)<2);% == min(abs(angle)));
                        loc = loc1;
%                         loc2 = find(abs(angle2)<0.1);% == min(abs(angle)));
%                         loc = [loc1;loc2];
                        if isempty(loc)% do I need to skip this point?Yes
%                             loc = find(abs(angle)== min(abs(angle)));
                            fprintf('here!\n');
                            continue;
                        end
                        dis = sqrt((x0-xr(loc)).^2+(y0-yr(loc)).^2+(z0-zr(loc)).^2);
                        thick(ind,1) = tmpbran(j);
                        if sum(dis<thickness(i,2)+11) == 0
                            thick(ind,2) = mean(dis);
                        else
                            thick(ind,2) = mean(dis((dis<thickness(i,2)+11)));
                        end
                        thick(ind,3) = j;
                        thick(ind,4)=i;
                        ind = ind+1;
                    end
                end
                if ~isempty(thick)
                    allthickness{index} = thick;
                    index = index+1;
                end
%             end
        end
end
toc

%%
save('allthickness.mat','allthickness')
% load allthickness
%% bran thickness
clear branthickness
clear branthicknessid
for i =1:length(allthickness)
    thick = allthickness{i}; 
    branthicknessid(i,1) = thick(1,4);
    branthickness(i,1) = min(thick(:,2))*0.65;
end
%% plot the branches based on the thickness center part
clf
color1=jet(11); % set colormap
 %%
cor = [];
nbins = linspace(1,10,12);
for i=1:3
    group=(nbins(i)<branthickness).*(branthickness<nbins(i+1));
    group=find(group);
    temp=[];
    if ~isempty(group)
        for j=1:size(group,1)
            temp=[temp;bran{branthicknessid(group(j))}];
            cor = [cor;branthicknessid(group(j))];
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
% cmap = colormap(color1) ;
% cbh = colorbar ; %Create Colorbar
% cbh.TickLabels = num2cell(linspace(1,10,6));

FigureFormat
xlim([1 datas(3)])
ylim([1 datas(1)])
zlim([1 datas(2)])
% saveas(gcf,'branch thickness_thick','epsc');
% saveas(gcf,'branch thickness_thick','pdf');
% saveas(gcf,'branch thickness_thick','png');
%% Correlation: branch orientation vs branch thickness analysis_center

figure
%
clf
cor = [];
for i=1:11
    group=(nbins(i)<branthickness).*(branthickness<nbins(i+1));
    group=find(group);
    temp=[];
    if ~isempty(group)
        for j=1:length(group)    
            cor = branthicknessid(group(j));
            [theta,phi]=myangle(connect(cor,1),connect(cor,2));
            temp(j,:)=real([theta,phi]);
        end
        if ~isempty(temp)
            p=polarplot(deg2rad(temp(:,1)),abs((temp(:,2)*180/pi)),'o','MarkerSize',6);
            p.Color = color1(i,:);
            hold on
        end
    end    
 
end

a=gca;
 a.LineWidth=2;
 a.FontSize=12;
 set(gcf,'color','white')
cmap = colormap(color1) ;
cbh = colorbar ; %Create Colorbar
cbh.TickLabels = num2cell(linspace(1,10,11));
% saveas(gcf,'Center region branch thickness-orientation','epsc')
%  saveas(gcf,'Center region branch thickness-orientationr','pdf')
%  saveas(gcf,'Center region branch thickness-orientation','png')
%% Branch thickness histogram

figure
nbins = 0:1:10;
hlc= histofrequency(branthickness,nbins);
x=hlc(:,1)-0.5;
bar(x,hlc(:,2))
% plot(hlc(:,1),hlc(:,2),'LineWidth',1)
% x=hlc(:,1);
% bar(x,hlc(:,2))
colormap jet
set(gcf,'color','white')
xlabel('Branch thickness');ylabel('Frequency')
% xticks([0 20 40 60 80 100])
% yticks([0,0.05,0.1,0.15,0.2])
% ylim([0,0.23])
ylim([0,0.3])
xlim([0,10])
% saveas(gcf,'Histogram thickness','epsc')
%  saveas(gcf,'Histogram thickness','pdf')
%  saveas(gcf,'Histogram thickness','png')
%%
% clf
% for i =1:length(allthickness)
%     thick = allthickness{i}; 
%     xx = thick(:,3);
% %     a = find(abs(thick(:,2))== min(thick(:,2)),1);
%     a = round(size(thick,1)/2);
%     xx = xx - thick(a,3);
%     
%     t = thick(:,2);
%     if (branlen(i)>40)&&(branlen(i)<10)
%         continue
%     else
% 
%     plot(xx,t)
%     hold on
%     end
%     if (max(xx)>40)
%         i
%     end
% end


%% fit all branch profiles

allfit = [];
clf
index = 1;
branlenTemp = [];
branlenThick = [];
OutLen = [];
OutThickness = [];
for i =1:size(allthickness,2)
    thick = allthickness{i};
    t = thick(:,2)*0.65;
    numbran = thick(1,4);
    allbran = bran{numbran};
    connectBran = connect(numbran,:);
    tmpbran = allbran;%allbran(thick(:,3));
    dis = 0;
    brancor = [];
    brancor(1) = 0;
    for j =1:length(tmpbran)-1
        [x1,y1,z1] = ind2sub(datas,tmpbran(j));
        [x2,y2,z2] = ind2sub(datas,tmpbran(j+1));
        dis = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2)*0.65+dis;
        brancor(j+1) = dis;
    end
    if sum(isnan(t))>0
        continue
    else
    if length(t)>1%original set
       
        a = find(abs(thick(:,2))== min(thick(:,2)),1);
        z = (brancor-brancor(a))';
        if (max(abs(z))>25) || ((dis)>35) || ((dis)<2)
            i
            continue
        else
        startDis = [vecdist(connectBran(1),tmpbran(1)),vecdist(connectBran(1),tmpbran(end))];
        EndDis = [vecdist(connectBran(2),tmpbran(1)),vecdist(connectBran(2),tmpbran(end))];
% [xs, ~] = sort(startDis);
        sDis = min(startDis);
        eDis = min(EndDis);
        branlenTemp(index) = dis+sDis*0.65+eDis*0.65;
       
        branlenThick(index) = branthickness(i);
        OutLen = [OutLen;z/dis];
        OutThickness = [OutThickness;t/branthickness(i)];
        index = index+1;
%         ft2 = fittype({'1','x^2'});
%         p = fit(z,t,ft2);
%         y1 = p(z);
%         plot(z,t,'.')
%         hold on
%         plot(z,y1)
%         allfit = [allfit;p.b];
%         legend(gca,'off')
%          hold on
%          ylim([0,15])
        end
    end
    end
end

%%
saveas(gcf,'profile','epsc')
 saveas(gcf,'profile','pdf')
 saveas(gcf,'profile','png')
%%
excelLenThick = [OutLen,OutThickness];
plot(OutLen,OutThickness,'.')
xlswrite('NormalizedLenThick',excelLenThick)
%%
excelLenThickPoint =  [branlenTemp;branlenThick]';
xlswrite('LenvsThick',excelLenThickPoint)
%%
clf
x = branlenTemp';
X = [ones(length(x),1) x];
y = branlenThick';
b = X\y;
xrange = min(x):0.1:max(x);
Xrange = [ones(length(xrange),1) xrange'];
yCalc2 = Xrange*b;
scatter(branlenTemp,branlenThick,'.','b','LineWidth',1)
hold on
plot(xrange,yCalc2,'--','LineWidth',1)
% grid on
ylim([0,10])

xlabel('Branch thickness (\mu m)')
ylabel('Branch length (\mu m)')
% saveas(gcf,'branlenvsThickness','epsc')
%  saveas(gcf,'branlenvsThickness','pdf')
%  saveas(gcf,'branlenvsThickness','png')
hold on
load allChainlength
load allChainthickness
allChainthickness(allChainlength>50) = [];
allChainlength(allChainlength>50) = [];
x = allChainlength;
X = [ones(length(x),1) x];
y = allChainthickness;
b = X\y;
xrange = min(x):0.1:max(x);
Xrange = [ones(length(xrange),1) xrange'];
yCalc2 = Xrange*b;
scatter(allChainlength,allChainthickness,'r','LineWidth',0.5)
hold on
plot(xrange,yCalc2,'--','LineWidth',1)
legend('All branches','Fitted curve','Chain','fitted','Location','best');