function [connect,newnode,bran]=merge_n(c,n,b)
% merge connected nodes 
global datas
connect=c;
bran=b;
node=n;
change=0;
newnode=node(1,1);
for i=2:size(node,1)
    nearnode=newnode(max(1,i-100):end);  %check nodes before this one
      n=neighbor(node(i,1),datas);
      num=sum(double(ismember(nearnode,n)));
      if  num==1
          cb = find(connect(:,1)==node(i,1));
          if isempty(cb)
          else
          connect(cb,1)=nearnode(ismember(nearnode,n));
          temp = reshape(bran{cb},[],1);
          bran{cb}=[node(i,1);temp];
          change=change+1;
          end
          cb = find(connect(:,2)==node(i,1));
          if isempty(cb)
          else
          connect(cb,2)=nearnode(ismember(nearnode,n));
          temp = reshape(bran{cb},[],1);
          bran{cb}=[temp;node(i,1)];
          change=change+1;
          end
      else
          newnode=[newnode,node(i,1)];
      end
  
end
newnode = newnode';
for i = 1:length(newnode)
    numcon = length(find(connect(:,1)==newnode(i)))+length(find(connect(:,2)==newnode(i)));
    newnode(i,2) = numcon;
end
end