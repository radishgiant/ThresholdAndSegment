function T=Cluster_Kittler(X,Trange,ClusterNum)
% reference:J. Kittler and J. Illingworth, "Minimum error thresholding," Pattern
% Recogn.19, 41¨C47,1986.
hg=hist(reshape(X,[],1),Trange);
hg=hg./sum(hg);
Ti=0;
JT=zeros(length(Trange),1);
for T=Trange
    Ti=Ti+1;
    P(1)=sum(hg(Trange<=T));
    P(2)=sum(hg(Trange>T));
    mean(1)=sum(hg(Trange<=T).*Trange(Trange<=T))/P(1);
    mean(2)=sum(hg(Trange>T).*Trange(Trange>T))/P(2);
    theta1=(Trange-mean(1)).^2.*hg;
    theta2=(Trange-mean(2)).^2.*hg;
    cov(1)=sum(theta1(Trange<=T))/P(1);
    cov(2)=sum(theta2(Trange>T))/P(2);
    JT(Ti)=1+2*(P(1)*log10(sqrt(cov(1)))+P(2)*log10(sqrt(cov(2))))-2*(P(1)*log10(P(1))+P(2)*log10(P(2)));
end
[~,Tind]=min(JT);
T=Trange(Tind);
end