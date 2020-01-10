function newbran = branch_order(bran,connect,datas)
newbran = cell(1);
for i = 1:length(bran)
    newTempbran = [];
    tmpbran= bran{i};
    startnode = connect(i,1);
    n = neighbor(startnode,datas);
    nextNode = intersect(n,tmpbran);
    if isempty(nextNode)
       clear dist
       for j = 1:length(tmpbran)
           dist(j) = vecdist(startnode,tmpbran(j));
       end
       ind = find(dist == min(dist));
       nextNode = tmpbran(ind);
    end
    newTempbran = [newTempbran,nextNode];
    flag = 1;
    while flag<length(tmpbran)
        flag = flag + 1;
        n = neighbor(nextNode,datas);
        curNode = setdiff(intersect(n,tmpbran),newTempbran);
        if ~isempty(curNode)
            newTempbran = [newTempbran,curNode];
            nextNode = curNode;
        else
            restBran = setdiff(tmpbran,newTempbran);
            clear dist
            for j = 1:length(restBran)
                dist(j) = vecdist(nextNode,restBran(j));
            end
            ind = find(dist == min(dist));
            newTempbran = [newTempbran,restBran(ind)];
            nextNode = restBran(ind); 
        end
    end
    newbran{i} = newTempbran';
end