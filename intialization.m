function [P,mean0,cov0,alpha]=intialization(Xd)
% intialization prior distrubtion,mean ,cov,and parameter alpha  that defines the range
% around in which pixels cannot be easily identified as either changed or unchanged.
rng('shuffle')
sizex=length(Xd);
alpha=randperm(9,1)/10;
Md=(max(Xd,[],1)-min(Xd,[],1))./2;
Tc=Md.*(1+alpha);
Tn=Md.*(1-alpha);
Sc=Xd(Xd>Tc);
Sn=Xd(Xd<Tn);
Pc=length(Sc)/sizex;
Pn=length(Sn)/sizex;
P=[Pc,Pn];
mean0=[mean(Sc),mean(Sn)];
cov0=[var(Sc),var(Sn)];
end