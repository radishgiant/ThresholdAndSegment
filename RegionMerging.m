function Rlabel=RegionMerging(varargin)
% post-processing over segment Image based RAG region merging
% input:
% I by M x N gray Image
% label: initial segment result
% Treshold: interation temination condion default is 1000
% RegionNum: the Region Number after Mereging defalut is 2(forground and
% background
% output:
% Rlable: merging result
% reference:
% writed by radishgiant


[I,Rlabel,Treshold,RegionNumber] = parse_inputs(varargin{:});
minsim=Treshold-1;
K=length(unique(Rlabel));
[matrixSize(1),matrixSize(2),channel]=size(I);
X=reshape(I,matrixSize(1)*matrixSize(2),channel);
Vset=cell(K,1);
Pset=cell(K,1);
delnode=[];
labelstart=min(Rlabel(:));
if labelstart==0
    Rlabel=Rlabel+1;
elseif labelstart>1
    Rlabel=Rlabel-labelstart+1;
end
for i=1:K
    [r,c]=find(Rlabel==i);
    Vset{i}=[r,c];
    Pset{i}=X(sub2ind(matrixSize,r,c),:);
end

while length(unique(Rlabel))>RegionNumber || minsim(1)<=Treshold
Adj=imRAG(Rlabel);
sim=partionsimliar(Pset,Adj,delnode);
[minsim,simInd]=sort(sim,1,'ascend');
node1=Adj(simInd(1),1);
node2=Adj(simInd(1),2);
node=min(node1,node2);
nodem=max(node1,node2);
Rlabel(Rlabel==node1|Rlabel==node2)=node;
Pset{node}=[Pset{node1};Pset{node2}];
Pset{nodem}=[];
Vset{node}=[Vset{node1};Vset{node2}];
Vset{nodem}=[];
delnode=[delnode;nodem];
end




end

function [I,label,Treshold,RegionNumber] = parse_inputs(varargin)

narginchk(2,4);
I = varargin{1};
label=varargin{2};
validateattributes(I,{'numeric'}, {'2d' '3d'}, ...
    mfilename, 'I', 1);
validateattributes(label,{'numeric'}, {'2d'}, ...
    mfilename, 'label', 2);
if nargin<3
    Treshold=1000;
    RegionNumber=2;
elseif nargin<4
    Treshold=varargin{3};
    RegionNumber=2;
else
    Treshold=varargin{3};
    RegionNumber=varargin{4};
validateattributes(Treshold,{'numeric'}, {'scalar','real'}, ...
    mfilename, 'Treshold', 3);
validateattributes(RegionNumber,{'numeric'}, {'scalar','real','>=',2}, ...
    mfilename, 'RegionNumber', 4);
end

end