function SearchofIJ=SearchNodeConnection(I,J,PvI)
SortofIJ=[I,J;J,I];      %线路连接关系
SearchofIJ=[];
for i=1:length(PvI)  %搜索和发电机节点连接的线路组合
    SearchofIJ=[SearchofIJ;SortofIJ(SortofIJ(:,1)==PvI(i),:)];    
end
SearchofIJ=unique(SearchofIJ,'rows');
end