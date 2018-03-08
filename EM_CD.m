function [Tp,varargout]=EM_CD(Xd,MaxIteration, logliklihoodThreshold);
% referecne:[25] L. Bruzzone and D. F. Prieto, ¡°Automatic analysis of the difference image
% for unsupervised change detection,¡± IEEE Trans. Geosci. Remote Sens.,
% vol. 38, no. 3, pp. 1170¨C1182, May 2000
% input
%      Xd is X2- X1 of size [sample,1]
% output
%      Tp is treshold based EM
sizex=length(Xd);
% intialization
fprintf('---------------------------------intialization paramater----------------------------\n')
[P0,mean0,cov0,alpha]=intialization(Xd);
fprintf(' the intialization paramater Pw=  [%f,%f],mean0=[%.3f,%.3f],var0=[%.3f,%.3f],alpha= %.3f\n',P0(1),P0(2),mean0(1),mean0(2),cov0(1),cov0(2),alpha)
mean_new=mean0;
cov_new=cov0;
P_new=P0;
%  changed class EM is wi=1
%  unchanged class EM is wj=2
wi=1;
wj=2;
cont=0;
logliklihoodThresOld = 0;
val=100;
fprintf('---------------------------------start EM Interation Algrithom----------------------------\n')
while(val>logliklihoodThreshold )
    cont=cont+1;
    
    %         Pwc=Gaussmodel(Xd,mean_new(wi),cov_new(wi));
    %         Pwn=Gaussmodel(Xd,mean_new(wj),cov_new(wj));
    Pwc=gaussmf(Xd,[sqrt(cov_new(wi)),mean_new(wi)]).*((2*pi*cov_new(wi))^(-1/2));
    Pwn=gaussmf(Xd,[sqrt(cov_new(wj)),mean_new(wj)]).*((2*pi*cov_new(wj))^(-1/2));
    % statistics Xd's distribution
    Pw=P_new(wi).*Pwc+P_new(wj).*Pwn;
    logliklihoodThresNew=sum(log(P_new(wi).*Pwc)+log(P_new(wj).*Pwn));
    val = abs((logliklihoodThresNew - logliklihoodThresOld)/logliklihoodThresOld);
    logliklihoodThresOld =logliklihoodThresNew;
    if cont >MaxIteration
        fprintf('up to MaxIteration')
        break;
    end
    %update mean and var
    mean_new(wi)=sum(P_new(wi).*Pwc./Pw.*Xd)/sum(P_new(wi).*Pwc./Pw);
    cov_new(wi)=sum(P_new(wi).*Pwc./Pw.*(bsxfun(@minus,Xd,mean_new(wi)).^2))/sum(P_new(wi).*Pwc./Pw);
    mean_new(wj)=sum(P_new(wj).*Pwn./Pw.*Xd)/sum(P_new(wj).*Pwn./Pw);
    cov_new(wj)=sum(P_new(wj).*Pwn./Pw.*(bsxfun(@minus,Xd,mean_new(wj)).^2))/sum(P_new(wj).*Pwn./Pw);
    % update prior ditributation Pc
    P_new(wi)=sum(P_new(wi).*Pwc./Pw)./sizex;
    P_new(wj)=sum(P_new(wj).*Pwn./Pw)./sizex;
    fprintf(' the %d interation changed and unchanged prior distributation of is %.3f and %f alterlatly\n',cont,P_new(wi),P_new(wj))
    fprintf(' the %d interation means of is [%f,%f]\n',cont,mean_new)
    fprintf(' the %d interation variance of is [%f,%f]\n',cont,cov_new)
   
    
end
cn=cov_new(2);
cc=cov_new(1);
un=mean_new(2);
uc=mean_new(1);
Pc=P_new(1);
Pn=P_new(2);
a=cn-cc;
b=2*(un*cc-uc*cn);
c=uc^2*cn-un^2*cc-2*cn*cc*log(sqrt(cn)*Pc/sqrt(cc)/Pn);
Tp=roots([a,b,c]);
% Tp=double(solve((x-train.mean(1))^2/2/train.cov(1)-(x-train.mean(2))^2/2/train.cov(2)-log(sqrt(train.cov(2))*train.P(1)/sqrt(train.cov(1))*train.P(2))))
fprintf('cc=%f,cn=%f,uc=%f,un=%f\n',cc,cn,uc,un);
fprintf('a=%f,b=%f,c=%f',a,b,c);

if nargout>1
    train.cov=cov_new;
    train.mean=mean_new;
    train.P=P_new;
    train.Pw=Pw;
    varargout{1}=train;
end
end

