clf

PlotScatter(tempRingDots(i,:)',10,'r')
% FigureFormat
% hold on
% [y,x,z] = ind2sub(datas,ringDots);
% quiver3(round(mean(z)),round(mean(x)),round(mean(y)),n(3),n(1),n(2),30,'Color','k','LineWidth',1,'MaxHeadSize',0.2)
% 
hold on 
PlotScatter(re,10,'g')
% FigureFormat
% hold on
% [y,x,z] = ind2sub(datas,ringDots);
% quiver3(round(mean(x)),round(mean(y)),round(mean(z)),n(1),n(2),n(3),30,'Color','k','LineWidth',1,'MaxHeadSize',0.2)
%%
PlotScatter(trimskel,1,'b',0)
FigureFormat
hold on
quiver3(0,0,0,0,0,1,30,'Color','r','LineWidth',1,'MaxHeadSize',0.5)
%% Detect directions
figure
% test = ring7{2,1};
[Fit,errors] = planeFitting_v2(branTemp);
[xt,yt,zt] = ind2sub(datas,branTemp); % x,y,z instead of y,x,z here
PlotScatter(branTemp,10,'r',0)
% hold on
% PlotScatter(trimskel,1,'b',0)
hold on
quiver3(mean(zt),mean(yt),mean(xt),Fit(3),Fit(2),Fit(1),30,'Color','r','LineWidth',1,'MaxHeadSize',0.5)
hold on
quiver3(mean(zt),mean(yt),mean(xt),1,0,0,30,'Color','g','LineWidth',1,'MaxHeadSize',0.5)
hold on
quiver3(mean(zt),mean(yt),mean(xt),0,1,0,30,'Color','b','LineWidth',1,'MaxHeadSize',0.5)
hold on
quiver3(mean(zt),mean(yt),mean(xt),0,0,1,30,'Color','k','LineWidth',1,'MaxHeadSize',0.5)

FigureFormat
temp = vectorAngle(Fit)
%% 
ring7 = ring{4};
dis = [];
allAngle7 = [];
allRingSize7 = [];
loc7 = [];
for j =1:size(ring7,1)
    index = 0;
    tempRingDots = ring7{j,1};
    tempRing = ring7{j,3};
    tempRingDir = ring7{j,2};
    if ~isempty(tempRing)
        angle = [];
        ringSize = [];
    for i =1:length(tempRing)
        tempDis = tempRing{i};
        if max(tempDis)>5
            continue
        else
            index = index+1;
            loc7 = [loc7;[j,i]];
        dis = [dis;tempDis];
        fit = tempRingDir(:,i);
        fit(4) = 1;
        tempAngle = vectorAngle([fit(1),fit(2),fit(4)]);
        angle = [angle;tempAngle];
        
        [y,x,z] = ind2sub(datas,tempRingDots(i,:)');
        z0 = -(fit(1)*x+fit(2)*y+fit(4))/fit(3);
        if length(find(z0<=0))>0
           ringdots = [x,y,z];
        else
           ringdots = [x,y,z0];
        end
        s = ringSize_7E(ringdots);
        ringSize = [ringSize;s];
        end
    end
    ring7{j,4} = angle;
    ring7{j,5} = ringSize;
    ring7{j,6} = index;
    allAngle7 = [allAngle7;angle];
    allRingSize7 = [allRingSize7;ringSize];
    end
end
