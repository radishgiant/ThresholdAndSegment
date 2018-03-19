function sim=partionsimliar(Pset,Adj,delnode)
for i=1:length(delnode)
Pset{delnode(i)}=inf;
end
% all region's mean
fN=@(x)(size(x{:},1));
Nset=arrayfun(fN,Pset,'UniformOutput',0);
fmean=@(x)(sum(x{:},1)./size(x{:},1));
u=arrayfun(fmean,Pset,'UniformOutput',0);
node1=Adj(:,1);
node2=Adj(:,2);
u=cell2mat(u);
Nset=cell2mat(Nset);
sim=(u(node1,:)-u(node2,:)).^2.*Nset(node1).*Nset(node2)./(Nset(node1)+Nset(node2));
% sim(delnode)=inf(length(delnode),1);
end