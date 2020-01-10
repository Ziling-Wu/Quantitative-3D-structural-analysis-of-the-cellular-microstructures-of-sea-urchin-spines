function nextnode = connected_node(current,connect)
% find connect nodes in connect matrix   
a1 = find(current == connect(:,1));
nextnode = [];
if ~isempty(a1)
        for j = 1:length(a1)
        nextnode = [nextnode;connect(a1(j),2)];
        end
end
a2 = find(current == connect(:,2));
if ~isempty(a2)
    for j = 1:length(a2)
        nextnode = [nextnode;connect(a2(j),1)];
    end
end