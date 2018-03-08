function T=Entropy_Kapur(X,Trange);
hg=hist(X,Trange);
pg=hg./sum(hg);
Ti=0;
TrangeSpan=Trange(pg>0);
JT=zeros(length(TrangeSpan),1);
for Tk=TrangeSpan
    Ti=Ti+1;
    PT=sum(pg(Trange<=Tk));
    H=-pg./PT.*log10(pg./PT);
    isnanind=~isnan(H);
    Hf=sum(H(Trange<=Tk&isnanind));
    Hb=sum(H(Trange>Tk&isnanind));
    JT(Ti)=Hf+Hb;
end
[~,Tind]=max(JT);
T=TrangeSpan(Tind);
end