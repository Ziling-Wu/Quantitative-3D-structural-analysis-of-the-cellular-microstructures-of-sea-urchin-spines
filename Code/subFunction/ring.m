function ring(nodef,branf,connectf,np)
% np is interested node
global datas
clf
nc = nodef;
bc=branf;
cc = connectf;
branc=[];
for i=1:size(bc,2)
     branc=[branc;bc{i}];
end
PlotScatter(branc,3,[0.5,0.5,0.5],0)
hold on

for ii=1:length(np)
    n0=nc(np(ii),1) ;
    clear lp5 lp6 lp7
    [lp5,lp6,lp7,nb1,nb2]=findloop(cc,n0);
    % Plot the loop
    lp=lp7;

    cMap=hsv(size(lp7,1)+size(lp6,1));
    for j=3:size(lp,1)-1
    bran=[];
    for i=1:size(lp,2)-1
         bran=[bran;bc{findbran(lp(j,i),lp(j,i+1),cc)}];
    end
    bran=[bran;bc{findbran(lp(j,1),lp(j,size(lp,2)),cc)}];

    PlotScatter(bran,30,cMap(6,:),0)
    hold on

    PlotScatter(lp(j,:)',100,'g',0)
    hold on
    end

    lp=lp6;

    for j=1:size(lp,1)-1
    bran=[];
    for i=1:size(lp,2)-1
         bran=[bran;bc{findbran(lp(j,i),lp(j,i+1),cc)}];
    end
    bran=[bran;bc{findbran(lp(j,1),lp(j,size(lp,2)),cc)}];

%     PlotScatter(bran,30,cMap(5,:),0)

%     PlotScatter(lp(j,:)',50,'r',0)
    end


    hold on
    lp=lp5;

    for j=1:size(lp,1)
        bran=[];
        for i=1:size(lp,2)-1
             bran=[bran;bc{findbran(lp(j,i),lp(j,i+1),cc)}];
        end
        bran=[bran;bc{findbran(lp(j,1),lp(j,size(lp,2)),cc)}];

%         PlotScatter(bran,30,cMap(4,:),0)
%         PlotScatter(lp(j,:)',50,'r',0)
    end 
    [x,y,z]=ind2sub(datas,n0);
%     a = 50;
%     axis([y-a y+a x-a x+a z-a z+a])
    FigureFormat
end
a = 100;
    axis([z-a z+a y-a y+a x-a x+a])