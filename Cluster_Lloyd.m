function T=Cluster_Lloyd(X,Trange)
covX=var(X);
hg=hist(X,Trange);
pg=hg./sum(hg);
Ti=0;
TrangeSpan=Trange(pg>0);
JT=zeros(length(TrangeSpan),1);
for Tk=TrangeSpan
    Ti=Ti+1;
    P(1)=sum(pg(Trange<=Tk));
    P(2)=sum(pg(Trange>Tk));
    mean(1)=sum(pg(Trange<=Tk).*Trange(Trange<=Tk));
    mean(2)=sum(pg(Trange>Tk).*Trange(Trange>Tk));
    JT(Ti)=(mean(1)+mean(2))/2+covX/(mean(1)-mean(2))*log10(P(2)/P(1));
end
[~,Tind]=min(JT);
T=Trange(Tind);
end