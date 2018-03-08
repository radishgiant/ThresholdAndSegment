function T=Entropy_Yen(X,Trange);
hg=hist(X,Trange);
pg=hg./sum(hg);
Ti=0;
TrangeSpan=Trange(pg>0);
JT=zeros(length(TrangeSpan),1);
for Tk=TrangeSpan
    Ti=Ti+1;
    PT=sum(pg(Trange<=Tk));
    Cf=-log10(sum((pg(Trange<=Tk)./PT).^2));
    Cb=-log10(sum((pg(Trange>Tk)./(1-PT)).^2));
    JT(Ti)=Cf+Cb;
end
JT(isinf(JT))=0;
JT(isnan(JT))=0;
[~,Tind]=max(JT);
T=TrangeSpan(Tind);
end