clearvars
%% Read in original image and do segmentation

%% load images
startpage = 1;
endpage = 100;
style = '.tif';
folder = 'D:\Ting data\Seaurchinspine\20190302binary';
% folder = '.\100x100'
img = myreadin(startpage,endpage,style,folder);

% Convert img into array with 0 and 1
%{
for i = 1:size(img,1)
    
    for j = 1:size(img,2)
        
        for k = 1:size(img,3)
            
            if img(i,j,k) ~= 0
                
                img(i,j,k) = 1;
                
            end
        
        end
        
    end    
    
end    
%}

%%  Skeleton
tic
skel_c = Skeleton3D(img);
disp(toc);
%%
datas = size(skel_c);
global datas

PlotScatter(skel_c,1);axis on;
xlabel('x')
ylabel('y')
zlabel('z')

% center: x growth
% 
%% Network cleaning by trimming branch
tic
trimnode = node_identification(skel_c);
trimskel = branch_trim(trimnode, skel_c);
trimskel = Skeleton3D(trimskel);
trimnode = node_identification(trimskel);
time = toc;
disp(time)
%% Label: node and branch ,
tic
[bran,connect, node] = branch_sort2(trimskel,trimnode);
time = toc;
disp(toc)

%%

figure
for i = 1:length(bran)
corbran = bran{i};
PlotScatter(corbran,1,'r') ;
hold on
end
PlotScatter(node(:,1),10,'g') ;
a = gca;
gui_FigureFormat(a,120,20)
hold off

%% Node merging
% merge connected nodes
tic
[c2,n2,b2] = merge_connectednodes(connect,node,bran);
time = toc;
disp(time)
% delete branches
tic
clear c3 n3 b3
[c3,n3,b3] = merge_double(c2,n2,b2,5);
time = toc;
disp(time)
%
% merge nearby nodes
tic
[c4,n4,b4] = merge_nearby(c3,n3,b3,10,1);
time = toc;
disp(time)
% delete branches again
clear c5 n5 b5
tic
[c5,n5,b5] = merge_double(c4,n4,b4);
time = toc;
disp(time)
%
nodef = [];
for i = 1:size(n5,1)
    if n5(i,2)>2
        nodef = [nodef;n5(i,:)];
    end
end
%
connectf = c5;
branf = b5;
%%
figure
for i = 1:length(b5)
corbran = b5{i};
PlotScatter(corbran,1,'r') ;
hold on
end
PlotScatter(n5(:,1),10,'g') ;
a = gca;
gui_FigureFormat(a,120,20)
title(a,'Final Plot')
%% 
%PlotScatter(trimskel)
%FigureFormat
%% Analysis: Node, branch length, branch thickness, branch orientation,ring, chain and FFT of node

%% Node
% node type
hlc= histofrequency(nodef(:,2),3:1:7);
x=hlc(:,1)-0.5;
bar(x,hlc(:,2))
colormap jet
set(gcf,'color','white')
%% Output node type along with their nodal position
filename = 'D:\Ting data\Sea urchin spine code\NodeTypeDistribution.csv';
file = fopen(filename,'w');
Boundary = 15;
for j = 1:length(nodef)
[y1,x1,z1] = ind2sub(datas,nodef(j,1)); % Fix this, ind2sub actually reverses the x and y coordinates

if x1>Boundary && x1<datas(2)-Boundary && y1>Boundary && y1<datas(1)-Boundary && z1>Boundary && z1<datas(3)-Boundary
    
fprintf(file,'%i',x1);
fprintf(file,'%s',',');
fprintf(file,'%i',y1);
fprintf(file,'%s',',');
fprintf(file,'%i',z1);
fprintf(file,'%s',',');
fprintf(file,'%i',nodef(j,2));
fprintf(file,'\n');

end

end
fclose(file);

%% Ring
np=[10,20,100,200,800,1000,2000,4000];
PlotSkel = 0; % 0 means do NOT plot the entire skeleton with rings; 1 means plotting the entire skeleton with the rings.
ring(nodef,branf,connectf,np,PlotSkel)

%% Chain
% Remember to choose starting nodes inside the Find_Chain function
threshold = 5; % tune this parameter
SearchAngle = 45; % In degree
[allchain45,loc_line45] = Find_Chain2(connectf, threshold, SearchAngle);
PlotSkel = 0; % 0 means do NOT plot the entire skeleton with rings; 1 means plotting the entire skeleton with the rings.
PlotChain(trimskel,connectf,branf,allchain45,loc_line45,PlotSkel);

%% FFT
node = nodef;
nodedistribution = zeros(datas);
nodedistribution(node(:,1))=1;
fftnode = abs((ifftshift(fftn(fftshift(nodedistribution)))));%*(0.65^3);
fftnode = fftnode - min(fftnode(:));
fftnode = fftnode./max(fftnode(:));
%%
figure
thr = 0.05;
cor = find(fftnode>thr);
[a,b,value] = find(fftnode>thr);
[cor1,cor2,cor3] = ind2sub(datas,cor); % find nonzero component
scatter3(cor2(:,1),cor1(:,1),cor3(:,1),10,value,'filled');
colormap gray
axis image
set(gcf,'color','white')

%% Branch central thickness
% thickness middle position of branch
[thickness,newbran_thickness,newconnect_thickness] = branch_thickness(bran,connect,node,img);
filename = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\PositionThickness.csv';
file = fopen(filename,'w');
thickPos = thickness(:,1);
Boundary = 15;
[y1,x1,z1] = ind2sub(datas,thickPos);
for j = 1:length(thickness)
    
if x1(j)>Boundary && x1(j)<datas(2)-Boundary && y1(j)>Boundary && y1(j)<datas(1)-Boundary && z1(j)>Boundary && z1(j)<datas(3)-Boundary
fprintf(file,'%i',x1(j));
fprintf(file,'%s',',');
fprintf(file,'%i',y1(j));
fprintf(file,'%s',',');
fprintf(file,'%i',z1(j));
fprintf(file,'%s',',');
fprintf(file,'%i',thickness(j,2));
fprintf(file,'\n');
end

end
fclose(file);

%% Branch length and distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Here begins the length and distance calculation and output %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
branlen = [];
brandis = [];
N = 0;
filename = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\Skeleton.csv';
file = fopen(filename,'w');
for i = 1: length(b5)

% Branch length
tmpbran= b5{i};
dis = 0;

for j =1:length(tmpbran)-1
        dis = vecdist(tmpbran(j),tmpbran(j+1))+dis;
end
startDis = [vecdist(c5(i,1),tmpbran(1)),vecdist(c5(i,1),tmpbran(end))];
EndDis = [vecdist(c5(i,2),tmpbran(1)),vecdist(c5(i,2),tmpbran(end))];
% [xs, ~] = sort(startDis);
sDis = min(startDis);
eDis = min(EndDis);
branlen(i) = dis+sDis+eDis;

% Branch distance
[y1,x1,z1] = ind2sub(datas,c5(i,1));
[y2,x2,z2] = ind2sub(datas,c5(i,2));
connectivity = [x1,y1,z1,x2,y2,z2];
brandis(i) = vecdist(c5(i,1),c5(i,2));

%{
for j =1:length(tmpbran)-1
    [x1,y1,z1] = ind2sub(datas,tmpbran(j));
    [x2,y2,z2] = ind2sub(datas,tmpbran(j+1));
    dis = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2)+dis;
end
startDis = [vecdist(c5(i,1),tmpbran(1)),vecdist(c5(i,1),tmpbran(end)),...
        vecdist(c5(i,2),tmpbran(1)),vecdist(c5(i,2),tmpbran(end))];
[xs, ~] = sort(startDis);
branlen(i) = dis+sum(xs(1:2));

% Branch distance
[y1,x1,z1] = ind2sub(datas,c5(i,1));
[y2,x2,z2] = ind2sub(datas,c5(i,2));
connectivity = [x1,y1,z1,x2,y2,z2];
N = N + 1;
brandis(i) = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2);
%}

% Output the connectivity to a file
for k = 1:length(connectivity)
    if k == length(connectivity)
       fprintf(file,'%s',int2str(connectivity(k)));
    else
       fprintf(file,'%s',strcat(int2str(connectivity(k)),','));
    end
end
fprintf(file,'\n');
end
fclose(file);

% disp(N)
figure
histogram(branlen(:))
xlabel(gca,'Branch length(voxels)')
ylabel(gca,'Number');set(gcf,'color','white')

branlengthtmp = branlen(:);
filename = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\BranchPathLengthVoxel.csv';
file = fopen(filename,'w');
for i = 1:length(branlengthtmp)
    fprintf(file,'%8.3f',branlengthtmp(i));
    fprintf(file,'\n');
end
fclose(file);

branlengthtmp = brandis(:);
filename = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\BranchDisVoxel.csv';
file = fopen(filename,'w');
for i = 1:length(branlengthtmp)
    fprintf(file,'%8.3f',branlengthtmp(i));
    fprintf(file,'\n');
end
fclose(file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Here begins the branch position and lendis output %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
branlen = [];
brandis = [];
BranLenLessThanDis = [];
NodeLenLessThanDis = [];
LenGreater100 = [];
NodeLenGreater100 = [];
N = 0;
Boundary = 20;
index = 0;
filename_Conn = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\ConnBranLenBranDisVoxel.csv';
file_Conn = fopen(filename_Conn,'w');
filename_Bran = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\BranBranLenBranDisVoxel.csv';
file_Bran = fopen(filename_Bran,'w');
for i = 1: length(b5)
tmpbran= b5{i};
dis = 0;
[y1,x1,z1] = ind2sub(datas,c5(i,1));
[y2,x2,z2] = ind2sub(datas,c5(i,2));

if (x1>Boundary && x1<datas(2)-Boundary && y1>Boundary && y1<datas(1)-Boundary && z1>Boundary && z1<datas(3)-Boundary) || ...
        (x2>Boundary && x2<datas(2)-Boundary && y2>Boundary && y2<datas(1)-Boundary && z2>Boundary && z2<datas(3)-Boundary)

connectivity = [x1,y1,z1,x2,y2,z2];
N = N + 1;
brandis(i) = vecdist(c5(i,1),c5(i,2));
for k = 1:length(connectivity)
    fprintf(file_Conn,'%s',strcat(int2str(connectivity(k)),','));
end
    
for j =1:length(tmpbran)-1
    dis = vecdist(tmpbran(j),tmpbran(j+1))+dis;
    
    [y1,x1,z1] = ind2sub(datas,tmpbran(j));
    [y2,x2,z2] = ind2sub(datas,tmpbran(j+1));
    BranConn = [x1,y1,z1,x2,y2,z2];
    for k = 1:length(BranConn)
        fprintf(file_Bran,'%s',strcat(int2str(BranConn(k)),','));
    end
        
end
% startDis = [vecdist(c5(i,1),tmpbran(1)),vecdist(c5(i,1),tmpbran(end)),...
%     vecdist(c5(i,2),tmpbran(1)),vecdist(c5(i,2),tmpbran(end))];
% [xs, ~] = sort(startDis);
startDis = [vecdist(c5(i,1),tmpbran(1)),vecdist(c5(i,1),tmpbran(end))];
EndDis = [vecdist(c5(i,2),tmpbran(1)),vecdist(c5(i,2),tmpbran(end))];
% [xs, ~] = sort(startDis);
sDis = min(startDis);
eDis = min(EndDis);
branlen(i) = dis+sDis+eDis;
% branlen(i) = dis+sum(xs(1:2));
fprintf(file_Bran,'%8.3f',branlen(i));
fprintf(file_Bran,'\n');

%{
if branlen(i) < brandis(i)
    index = index + 1;
    NodeLenLessThanDis(index,1) = c5(i,1);
    NodeLenLessThanDis(index,2) = c5(i,2);
    for j = 1:length(tmpbran)
        BranLenLessThanDis(index,j) = tmpbran(j);
    end
end
%}

if branlen(i) > 100
    index = index + 1;
    NodeLenGreater100(index,1) = c5(i,1);
    NodeLenGreater100(index,2) = c5(i,2);
    for j = 1:length(tmpbran)
        LenGreater100(index,j) = tmpbran(j);
    end
end
    

distant = brandis(i);
fprintf(file_Conn,'%8.3f',branlen(i));
fprintf(file_Conn,'%s',',');
fprintf(file_Conn,'%8.3f',distant);
fprintf(file_Conn,'\n');
end

end
fclose(file_Conn);
fclose(file_Bran);

%%
figure
NodeLenGreater100Unique = unique(NodeLenGreater100);
for i = 1:length(LenGreater100)
branNonZero = find(LenGreater100(i,:));
for j = 1:length(branNonZero)
corbran(j) = LenGreater100(i,branNonZero(j));
end
PlotScatter(corbran,10,'r') ;
hold on
end
PlotScatter(NodeLenGreater100Unique,10,'g');
hold off
%PlotScatter(node(:,1),10,'g') ;
%a = gca;
%gui_FigureFormat(a,120,20)
%hold off
%{
figure
NodeLenLessThanDisUnique = unique(NodeLenLessThanDis);
for i = 1:length(BranLenLessThanDis)
branNonZero = find(BranLenLessThanDis(i,:));
for j = 1:length(branNonZero)
corbran(j) = BranLenLessThanDis(i,branNonZero(j));
end
PlotScatter(corbran,10,'r') ;
hold on
end
PlotScatter(NodeLenLessThanDisUnique,10,'g');
%PlotScatter(node(:,1),10,'g') ;
a = gca;
gui_FigureFormat(a,120,20)
hold off
%}

%% Branch profile
bran = branf;
connect = connectf;
branlen = [];
brandis = [];
for i =1: length(bran)
    tmpbran= bran{i};
    dis = 0;
    for j =1:length(tmpbran)-1
        dis = vecdist(tmpbran(j),tmpbran(j+1))+dis;
    end
    % startDis = [vecdist(connect(i,1),tmpbran(1)),vecdist(connect(i,1),tmpbran(end)),...
    %     vecdist(connect(i,2),tmpbran(1)),vecdist(connect(i,2),tmpbran(end))];
    % [xs, ~] = sort(startDis);
    % branlen(i) = dis+sum(xs(1:2));
    startDis = [vecdist(c5(i,1),tmpbran(1)),vecdist(c5(i,1),tmpbran(end))];
    EndDis = [vecdist(c5(i,2),tmpbran(1)),vecdist(c5(i,2),tmpbran(end))];
    sDis = min(startDis);
    eDis = min(EndDis);
    branlen(i) = dis+sDis+eDis;
    brandis(i) = vecdist(connect(i,1),connect(i,2));
end

% histogram(branlen)
%
ratio =[];
index = 1;
for i =1:length(brandis)
%     if branlen(i)>10
%         if branlen(i)<30
            ratio(index) = branlen(i)/brandis(i);
            index = index+1;
%         end
%     end
end
%
histogram(ratio)
xlim([0.8,2])

%%% calculate branch profile
%%%
%%% This branch profile code uses sphere to fit the local branch thickness
%%% Minimal sphere fitting is applied? 
%%% (Ref: Statistical Analysis of the Local Strut Thickness of Open Cell Foams)
%%% 
cont = edge3(img,'approxcanny',0.6);
cont = permute(cont,[2 1 3]);
tic
% [y,x,z]=meshgrid(1:datas(2),1:datas(1),1:datas(3));
thickness = cell(1);
Coordx = cell(1);
Coordy = cell(1);
Coordz = cell(1);
index = 1;
conn=connectf;
% branForThick=cell2mat(branf);
width = 15; %this means we only care about the branches that are not closed to the margin of the model
Lundesired = 0;
Rdesired = 0;
for i = 1:size(conn,1)
    
    tmpbran = cell2mat(branf(i));
    [y1,x1,z1] = ind2sub(datas,conn(i,1));
    [y2,x2,z2] = ind2sub(datas,conn(i,2));
    
    l = length(tmpbran);
    
    if l<10 || l>100  % Not too short and not too long
        Lundesired = Lundesired + 1;
        continue
    elseif (ratio(i)>=1)&&(ratio(i)<=2.5) %not too curve branch, change this parameter everytime
        Rdesired = Rdesired + 1;
        [yy,xx,zz] = ind2sub(datas,tmpbran);
        minx = min(xx);miny = min(yy);minz = min(zz);
        maxx = max(xx);maxy = max(yy);maxz = max(zz); % Find minimum and maximum: the bounding box of this branch
        flag = minx>width && maxx<datas(2)-width &&miny>width && maxy<datas(1)-width &&minz>width && maxz<datas(3)-width;
        
        if flag
            
            ind = 1;
            thick = [];
            
            for j = 1:length(tmpbran)
                
                w = 0;
                temp = 0;
                [y0,x0,z0] = ind2sub(datas,tmpbran(j));
                [y,x,z]=meshgrid((y0-15):(y0+15),(x0-15):(x0+15),(z0-15):(z0+15));
                
                while temp<=10 && w<15
                    % Find the just right size of bounding box (10 intersecting points, this parameter can be tuned as well, 1 is also worth a try)
                    
                    w = w+1; % w will finally become the thickness of the current position
                    % Sphere = (((x-x0).^2+(y-y0).^2+(z-z0).^2) <= w^2);
                    region = cont((x0-15):(x0+15),(y0-15):(y0+15),(z0-15):(z0+15)).*(((x-x0).^2+(y-y0).^2+(z-z0).^2) <= w^2);
                    % region = cont(x0-w:x0+w,y0-w:y0+w,z0-w:z0+w);
                    % region = cont(x0-w:x0+w,y0-w:y0+w,z0-w:z0+w);
                    % Expanding the search region (A retangle region, different from the spherical region in the main function)
                    temp = length(find(region == 1));
                    
                end    
                if w == 14
                    fprintf('%d',w);
                end    
                % Output the thickness and the position information
                thick(ind,1) = tmpbran(j); % Index of coordinate of the current location
                thick(ind,2) = w; % The thickness at the current location
                thick(ind,3) = j; % The location of the current position of this branch 
                thick(ind,4) = i; % Indicate which branch
                ind = ind+1;
            end
            
            thickness{index} = thick;
            Coordx{index} = (x1+x2)/2;
            Coordy{index} = (y1+y2)/2;
            Coordz{index} = (z1+z2)/2;
                
            index = index+1;
            
        end
    end
end
toc
%{
for i = 1:size(c3,1)
    tmpbran = b3{i};
    [y1,x1,z1] = ind2sub(datas,c3(i,1));
    [y2,x2,z2] = ind2sub(datas,c3(i,2));
    r1 = [x2-x1,y2-y1,z2-z1]; % This is the node-to-node vector
    width = 15; %this means we only care about the branches that are not closed to the margin of the model
    l = length(tmpbran);
    if l<10 || l>50 % Not too short and too long
        continue
    else
        if (ratio(i)>=1)&&(ratio(i)<=1.2) %not too curve branch
            [xx,yy,zz]= ind2sub(datas,tmpbran); % Get all the coordinates of one branch
            minx = min(xx);miny = min(yy);minz = min(zz);
            maxx = max(xx);maxy = max(yy);maxz = max(zz); % Find minimum and maximum: the bounding box of this branch
            flag = minx>width && maxx<datas(1)-width &&miny>width && maxy<datas(2)-width &&minz>width && maxz<datas(3)-width;
            if flag
                ind = 1;
                thick=[];
                for j = 1:length(tmpbran)
                    [y0,x0,z0] = ind2sub(datas,tmpbran(j));
                    % tempsphere = cont.*(((x-y0).^2+(y-x0).^2+(z-z0).^2) <= 20^2); 
                    tempsphere = cont.*(((x-x0).^2+(y-y0).^2+(z-z0).^2) <= 8^2); 
                    % What does the above line mean??? 
                    % This means, if the later part, or the sphere radius?, is less or equal than 20, 
                    % use the value for multiplicatyion, zero otherwise.
                    % The later part means: The distance between any points
                    % in the model and the point on the skeletal branch
                    % less than 20 voxels. The cont array helps checking if
                    % the point on the model.
                    
                    cor = find(tempsphere); % See how many are not 
                    if isempty(cor)
                        continue
                    else
                        [yr,xr,zr] = ind2sub(datas,cor);
                        c = [xr-x0,yr-y0,zr-z0];
                        clear angle1
                        angle1 = abs(acos(c*r1'/norm(r1)./vecnorm(c,2,2))*180/pi)-90; 
                        % Unit in degree
                        % Calculate the angle between the node-to-node vector and all other branch point-to-cor vectors
%                       angle2 = abs(acos(c*r2'/norm(r2)./vecnorm(c,2,2))*180/pi)-90;
                        loc1 = find(abs(angle1)<0.8);% == min(abs(angle)));
                        % This step finds the points that can form a vector
                        % with branch point and this vector is going to be
                        % nearly perpendicular to the node-to-node vector
                        loc = loc1;
%                         loc2 = find(abs(angle2)<0.1);% == min(abs(angle)));
%                         loc = [loc1;loc2];
                        if isempty(loc)% do I need to skip this point?Yes
%                             loc = find(abs(angle)== min(abs(angle)));
                            fprintf('here!\n');
                            continue;
                        end
                        thick(ind,1) = tmpbran(j);
                        thick(ind,2) = mean(sqrt((x0-xr(loc)).^2+(y0-yr(loc)).^2+(z0-zr(loc)).^2));
                        thick(ind,3) = j;
                        thick(ind,4) = i;
                        ind = ind+1;
                    end
                end
                if ~isempty(thick)
                    thickness{index} = thick;
                    Coordx{index} = (x1+x2)/2;
                    Coordy{index} = (y1+y2)/2;
                    Coordz{index} = (z1+z2)/2;
                    
                    index = index+1;
                    
                end
                end
            end
        end
        
end
%}

%
figure
filename = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\NearlyAllBranThickness.csv';
file = fopen(filename,'w');
CoordxArray = cell2mat(Coordx);
CoordyArray = cell2mat(Coordy);
CoordzArray = cell2mat(Coordz);
for i =1:length(thickness)
    
    thick = thickness{i};
    xx = thick(:,3);
    a = find(abs(thick(:,2))== min(thick(:,2)),1);
    xx = xx - a; % What does this mean??? This means move the original point to the minimum thickness position
    % xx is not the actual position or coordinates of the branch points,
    % but is voxel
    t = thick(:,2);
    
    % Output the thickness for each branch here
    % xx is the position array of the thicknesses
    fprintf(file,'%i',i);
    fprintf(file,'%s',',');
    for j = 1:length(xx)
    if j ~= length(xx)
    fprintf(file,'%8.3f',xx(j));
    fprintf(file,'%s',',');
    else
    fprintf(file,'%8.3f',xx(j));
    end
    end
    fprintf(file,'\n');
    
    xm = CoordxArray(i);
    fprintf(file,'%8.3f',xm);
    fprintf(file,'%s',',');
    ym = CoordyArray(i);
    fprintf(file,'%8.3f',ym);
    fprintf(file,'%s',',');
    zm = CoordzArray(i);
    fprintf(file,'%8.3f',zm);
    fprintf(file,'%s',',');
    
    for k = 1:length(t)
    if k ~= length(t)
    fprintf(file,'%8.3f',t(k));
    fprintf(file,'%s',',');
    else
    fprintf(file,'%8.3f',t(k));
    end
    end
    fprintf(file,'\n');
    
    plot(xx,t,'.')
    hold on
end
fclose(file);

% fit all branches 
%{
allfit = [];
clf
indx=1;
x1 = linspace(-40,40);
branlen  = [];
for i =1:size(thickness,2)
    thick = thickness{i};
    numbran = thick(1,4);
    allbran = b3{numbran};
    tmpbran = allbran(thick(:,3));
    dis = 0;
    brancor = [];
    brancor(1) = 0;
    for j =1:length(tmpbran)-1
        [x1,y1,z1] = ind2sub(datas,tmpbran(j));
        [x2,y2,z2] = ind2sub(datas,tmpbran(j+1));
        dis = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2)+dis;
        brancor(j+1) = dis;
    end
    branlen(index) = dis;
    index = index+1;
    t = thick(:,2);
    if length(t)>5%original set
%         a = round(max(thick(:,3))/2);
        a = find(abs(thick(:,2))== min(thick(:,2)),1);
        z = (brancor-brancor(a))';
%         xx = thick(:,3);
%         xx = xx - a;
%         z = xx*0.65;
        mint =find(abs(thick(:,2))== min(thick(:,2)),1);
        ft2 = fittype({'1','x^2'});
        p = fit(z,t,ft2);
        y1 = p(z);
%         p = polyfit(z,t,2);
%         y1 = polyval(p,z);
        plot(z,t,'.')
        hold on
        plot(z,y1)
        allfit = [allfit;p.b];
        legend(gca,'off')
         hold on
         ylim([2,20])
%          pause(1)
    end
end
set(gcf,'color','white')
%}
%% Branch orientation
% these lines only calcualte orientation of each branch
% in the paper, we care more about the correlation between branch
% orientation and branch length, branch thickness
clear connect
connect = connectf;
clf
% clear group temp
% color=jet(11);

% Modify or check the myangle function everytime using it!!!
filename = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\BranMidPosThetaPhixdir.csv';
file = fopen(filename,'w');
for i=1:size(connect,1)
    
    [y1,x1,z1] = ind2sub(datas,connect(i,1));
    [y2,x2,z2] = ind2sub(datas,connect(i,2));
    
    if (x1>Boundary && x1<datas(2)-Boundary && y1>Boundary && y1<datas(1)-Boundary && z1>Boundary && z1<datas(3)-Boundary) || ...
        (x2>Boundary && x2<datas(2)-Boundary && y2>Boundary && y2<datas(1)-Boundary && z2>Boundary && z2<datas(3)-Boundary)
    
    [theta,phi,z,coordm]=myangle(connect(i,1),connect(i,2));
    fprintf(file,'%8.3f',coordm(1));
    fprintf(file,'%s',',');
    fprintf(file,'%8.3f',coordm(2));
    fprintf(file,'%s',',');
    fprintf(file,'%8.3f',coordm(3));
    fprintf(file,'%s',',');
    fprintf(file,'%8.3f',theta);
    fprintf(file,'%s',',');
    fprintf(file,'%8.3f',phi);
    fprintf(file,'\n');
    % temp(i,:)=real([theta,phi]);
    
    end
    
end
fclose(file);
%
p=polarplot(temp(:,1)/pi*180,temp(:,2)/pi*180,'o');

%% Interbranch angles
nc = nodef;
cc = connectf;
filename3N = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\N3PosConnInterAngles.csv';
filename4N = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\N4PosConnInterAngles.csv';
file3N = fopen(filename3N,'w');
file4N = fopen(filename4N,'w');
Boundary = 15;
[ynode,xnode,znode] = ind2sub(datas,nc(:,1));
for i=1:size(nc,1)
    
    [y0,x0,z0] = ind2sub(datas,nc(i,1));
    
    if x0>Boundary && x0<datas(2)-Boundary && y0>Boundary && y0<datas(1)-Boundary && z0>Boundary && z0<datas(3)-Boundary
    
    temp = cc(cc(:,1)== nc(i,1),2); % Find the node index that is connected to the current node from the first column
    temp = [temp;cc(cc(:,2)==nc(i),1)]; % Find the node index that is connected to the current node from the second column
    for j=1:length(temp)
        angle_c(i,j)=interangle(nc(i,1),temp(j),temp(mod(j,length(temp))+1))*180/pi; % The angle between the i-th node and the j-th node
    end
    
    if nc(i,2) == 3
    
        fprintf(file3N,'%8.3f',xnode(i));
        fprintf(file3N,'%s',',');
        fprintf(file3N,'%8.3f',ynode(i));
        fprintf(file3N,'%s',',');
        fprintf(file3N,'%8.3f',znode(i));
        fprintf(file3N,'%s',',');
        fprintf(file3N,'%i',nc(i,2));
        fprintf(file3N,'%s',',');
    
        for k = 1:3
        
            fprintf(file3N,'%8.3f',angle_c(i,k));
            if k ~= 3
                fprintf(file3N,'%s',',');
            end
        
        end
    
        fprintf(file3N,'\n');
    
    elseif nc(i,2) == 4
        
        fprintf(file4N,'%8.3f',xnode(i));
        fprintf(file4N,'%s',',');
        fprintf(file4N,'%8.3f',ynode(i));
        fprintf(file4N,'%s',',');
        fprintf(file4N,'%8.3f',znode(i));
        fprintf(file4N,'%s',',');
        fprintf(file4N,'%i',nc(i,2));
        fprintf(file4N,'%s',',');
    
        for k = 1:4
        
            fprintf(file4N,'%8.3f',angle_c(i,k));
            if k ~= 4
                fprintf(file4N,'%s',',');
            end
        
        end
        
        fprintf(file4N,'\n');
    
    end
    
    end
        
end
interAngle = angle_c;
fclose(file3N);
fclose(file4N);

% 3N node
[x,hlcm,hlcs,hlcmi,hlcl] = Interbranchangle3N(interAngle);
XX = x(:);
MEAN = hlcm(:);
SMALL = hlcs(:); 
MEDIUM = hlcmi(:);
LARGE = hlcl(:);
filename = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\N3AngleMeanSML.csv';
file = fopen(filename,'w');
for i = 1:length(XX)
    fprintf(file,'%8.3f',XX(i));
    fprintf(file,'%s',',');
    fprintf(file,'%i',MEAN(i));
    fprintf(file,'%s',',');
    fprintf(file,'%i',SMALL(i));
    fprintf(file,'%s',',');
    fprintf(file,'%i',MEDIUM(i));
    fprintf(file,'%s',',');
    fprintf(file,'%i',LARGE(i));
    fprintf(file,'\n');
 end
fclose(file);
% 4N node
[x4,hlcm4,hlcs4,hlcl4] = Interbranchangle4N(interAngle);
XX = x4(:);
MEAN = hlcm4(:);
SMALL = hlcs4(:);
LARGE = hlcl4(:);
filename = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\N4AngleMeanSL.csv';
file = fopen(filename,'w');
for i = 1:length(XX)
    fprintf(file,'%8.3f',XX(i));
    fprintf(file,'%s',',');
    fprintf(file,'%i',MEAN(i));
    fprintf(file,'%s',',');
    fprintf(file,'%i',SMALL(i));
    fprintf(file,'%s',',');
    fprintf(file,'%i',LARGE(i));
    fprintf(file,'\n');
end
fclose(file);


% Planarity
% 
np = nodef;
cp = connectf;
filename = 'D:\Hongshun\Sea_Urchin\SUStructureAnalysisInvitedPaper\NewNewNewData_20190220\Growth8to10All\NoBoundaryData\PlanarityPosIndexV2.csv';
file = fopen(filename,'w');
Boundary = 15;
for i = 1:size(np)
    
    if np(i,2) == 3
        
        
        [y0,x0,z0] = ind2sub(datas,np(i,1));
    
        if x0>Boundary && x0<datas(2)-Boundary && y0>Boundary && y0<datas(1)-Boundary && z0>Boundary && z0<datas(3)-Boundary
        
        temp = cc(cc(:,1)== nc(i,1),2); % Find the node index that is connected to the current node from the first column
        temp = [temp;cc(cc(:,2)==nc(i),1)]; % Find the node index that is connected to the current node from the second column
        
        % Coords of the nodes that form a plane
        
        [ynode1,xnode1,znode1] = ind2sub(datas,temp(1));
        [ynode2,xnode2,znode2] = ind2sub(datas,temp(2));
        [ynode3,xnode3,znode3] = ind2sub(datas,temp(3));
        
        % Calculate the unit normal of the plane
        vnode = [x0,y0,z0];
        
        v1 = [xnode1-x0,ynode1-y0,znode1-z0];
        v1 = vnode + v1/norm(v1);
        
        v2 = [xnode2-x0,ynode2-y0,znode2-z0];
        v2 = vnode + v2/norm(v2);
        
        v3 = [xnode3-x0,ynode3-y0,znode3-z0];
        v3 = vnode + v3/norm(v3);
        
        v12 = [v2(1)-v1(1),v2(2)-v1(2),v2(3)-v1(3)];
        v23 = [v3(1)-v2(1),v3(2)-v2(2),v3(3)-v2(3)];
        
        vn = cross(v12,v23);
        vn = vn/norm(vn);
        
        vp = [x0-v1(1),y0-v1(2),z0-v1(3)];
        
        % Calculate planarity index
        D = abs(dot(vn,vp))/norm(vn);
        
        % Output
        fprintf(file,'%8.3f',x0);
        fprintf(file,'%s',',');
        fprintf(file,'%8.3f',y0);
        fprintf(file,'%s',',');
        fprintf(file,'%8.3f',z0);
        fprintf(file,'%s',',');
        fprintf(file,'%8.3f',D);
        fprintf(file,'\n');
        
        end
        
    end
        
end
fclose(file);