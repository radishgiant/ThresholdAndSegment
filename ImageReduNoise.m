function [X,varargout]=ImageReduNoise(I,varn,varargin)
global Ctable NeigborHoodSizeSet SigLevel1 SigLevel2 SigLevel3 NeigborHoodSizeInd NeigborHoodSizeSet3d
NeigborHoodSizeSet=[3,3;5,5;7,7;9,9;11,11;...
    13,13;15,15;17,17;19,19];
NeigborHoodSizeSet3d=[3,3,3;5,5,5;7,7,7];
NeigborHoodSizeInd=[1;2;3;4;5;6;7;...
    8;9;10;11;12];
SigLevel1=[0.7230;0.4566;0.3265;0.2577;0.2181;...
    0.1837;0.1787;0.1397;0.1247;0.4402;...
    0.2154;0.1282];
SigLevel2=[0.9483;0.5746;0.4286;0.3164;0.2654;...
    0.2228;0.1920;0.1688;0.1505;0.5527;...
    0.2609;0.1546];
SigLevel3=[1.2322;0.7192;0.5102;0.3868;0.3553;...
    0.2974;0.2557;0.2243;0.2000;0.6905;...
    0.3492;0.2051];
Ctable=table(NeigborHoodSizeInd,SigLevel1,SigLevel2,SigLevel3,...
    'VariableNames',{'NeigborHoodSize' 'SigLevel1' 'SigLevel2' 'SigLevel3'});
nargoutchk(1,2)
narginchk(2,6);
[NeigborHoodSize,SigLevel] = parse_inputs(varargin{:});
eval(['C=Ctable.',SigLevel,'(Ctable.NeigborHoodSize==NeigborHoodSize);'])

n=NeigborHoodSizeSet(NeigborHoodSize,1);
[M,N,C]=size(I);
X=zeros(M,N,C);
SNR=zeros(M*N,C);
tempI=zeros(M+(n-1),N+(n-1),C);
tempI((n-1)/2+1:M+(n-1)/2,(n-1)/2+1:N+(n-1)/2,:)=I;
I=tempI;clear tempI;
for rowi=(n-1)/2+1:M+(n-1)/2
    for coli=(n-1)/2+1:N+(n-1)/2
        rowstart=rowi-(n-1)/2;
        rowend=rowi+(n-1)/2;
        colstart=coli-(n-1)/2;
        colend=coli+(n-1)/2;
        for ci=1:C
            model=I(rowstart:rowend,colstart:colend,ci);
            Su=sum(model(:))./(n*n);
            S=sum((model(:)-Su).^2)./(n*n);
            if S<=(varn(ci)*(1+C))
                X(rowi-(n-1)/2,coli-(n-1)/2,ci)=Su;
            else
                % ¹À¼Æ²ÎÊý
                m1=Su;
                m2=sum(model(:).^2)./(n*n);
                m3=sum(model(:).^3)./(n*n);
                c1=m1;
                c2=m2-varn(ci);
                c3=m3-3*m1*varn(ci);
                beta=(c3-c1*c2)/(c2-c1^2);
                lamda=(c1*c3-c2*c2)/(c2-c1^2);
                u0=(beta-sqrt(beta^2-4*lamda))/2;
                u1=(beta+sqrt(beta^2-4*lamda))/2;
                P0=(u1-c1)/(u1-u0);
                P1=-(u0-c1)/(u0-u1);
                Tc=u0/2+u1/2+varn(ci)/(u1-u0)*log(P0/P1);
                SNR(sub2ind([M,N],rowi-(n-1)/2,coli-(n-1)/2),ci)=abs(u0-u1)/sqrt(varn(ci));
                % estimate X
                if I(rowi,coli)>Tc
                    X(rowi-(n-1)/2,coli-(n-1)/2,ci)=u1;
                else
                    X(rowi-(n-1)/2,coli-(n-1)/2,ci)=u0;
                end
            end
        end
    end
    
end
if nargout==2
    varargout{1}=SNR;
end
end
function [NeigborHoodSize,SigLevel] = parse_inputs(varargin)
global NeigborHoodSizeSet NeigborHoodSizeInd
pnames = { 'NeigborHoodSize','SigLevel'};
dflts =  {[13,13],'SigLevel1'};
[NeigborHoodSize,SigLevel] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});
validateattributes(NeigborHoodSize,{'numeric'}, {'real' 'row'});
SigLevelNames = {'SigLevel1','SigLevel2','SigLevel3'};
SigLevel = internal.stats.getParamVal(SigLevel,SigLevelNames,'''Method''');
if size(NeigborHoodSize,2)==2
    if~nnz(ismember(NeigborHoodSizeSet,NeigborHoodSize,'rows'))
        error('NeigborHoodSize is not legal')
    else
        NeigborHoodSize=NeigborHoodSizeInd(ismember(NeigborHoodSizeSet,NeigborHoodSize,'rows'));
    end
end
if size(NeigborHoodSize,2)==3
    if ~nnz(ismember(NeigborHoodSizeSet3d,NeigborHoodSize,'rows'))
        error('NeigborHoodSize is not legal')
    else
        Ind=find(ismember(NeigborHoodSizeSet3d,NeigborHoodSize,'rows'))+9;
        NeigborHoodSize=NeigborHoodSizeInd(Ind);
    end
end


end

