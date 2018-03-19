function [P_new,label]=Local_Yanowitz(I,varargin);
% local adaptive  treshold segment by Yanowitz
% input:
%       I is Image M by N
% writed by radishgiant
% github:https://github.com/radishgiant/ThresholdAndSegment.git
%reference:S. D. Yanowitz and A. M. Bruckstein, "A new method for image
% segmentation," Comput. Graph. Image Process. 46, 82¨C95 ,1989.
dbstop if error
if nargin<=1||isempty(varargin{1});
    hsize=[3,3];
end
if nargin<=2||isempty(varargin{2});
    MaxInterNum=15500;
end
if nargin<=3||isempty(varargin{3});
    InterTreshhold=10e-6;
end
if nargin<=4||isempty(varargin{4});
    GradTresh=20;
end
% I=double(I);
[M,N]=size(I);
% step1:smooth the image
h1=fspecial('average',hsize);
SI=imfilter(I,h1);
%step2: calculate gradiant map
[Fx,Fy]=gradient(double(SI));
F=sqrt(Fx.^2+Fy.^2);
%step3 :Laplacian Image
h2=fspecial('laplacian',0);
LI=imfilter(SI,h2);
%step4:sample the smoothed image at the places which the maximal
%gradiant-mask points
P_new=zeros(M,N);
P_new(LI==0)=I(LI==0);

% step5: interpolate the sampled gray level over the image
Residual=InterTreshhold+1;
InterNum=0;
while (Residual>InterTreshhold)
    if(InterNum>MaxInterNum)
        fprintf('up to MaxInterNum without diveregence');
        break;
    end
    InterNum=InterNum+1;
    P_last=P_new;
    R=imfilter(P_new,h2);
    P_new=P_new+R./4;
    Residual=mean(abs(P_new(:)-P_last(:)));
    
end

% step:6 segment the Image
bw=zeros(M,N);
bw(I>P_new)=255;%background
figure,imshow(bw);title('first segment result')
% step:7 validation progress
label=bwlabel(bw,4);
RGBLabel=label2rgb(label);
figure,imshow(RGBLabel);title('connected component');
lable_n=length(unique(label));
gradientmean=zeros(lable_n,1);
toglabel=zeros(lable_n,1);
for ci=0:lable_n-1
    temp=zeros(size(I));
    temp(label==ci)=255;
    eg=edge(temp);
    gradientmean(ci+1)=mean(F(eg==1));
    [egr,egc]=find(eg==1);
    [~,mingI]=min(F(eg>=1));% find the location of gradient of min value in eg
    mingr=egr(mingI);%find the location of gradient of min value in over image
    mingc=egc(mingI);
   
        nearborlabel=[mingr+1,mingc;mingr-1,mingc;mingr,mingc+1;mingr,mingc-1];
    nearborlogical=ones(4,1);
        if (mingr==1)
        nearborlogical(2)=0;
        end
        if (mingr==M)
        nearborlogical(1)=0;
        end
    if (mingc==1)
         nearborlogical(4)=0;
    end
    if mingc==N
       nearborlogical(3)=0;
    end
    nearborlabel=nearborlabel(nearborlogical==1,:);
    nearborlabel=label(sub2ind([M,N],nearborlabel(:,1),nearborlabel(:,2)));
    dlilabel=label(mingr,mingc);
    if nnz(nearborlabel~=dlilabel)
        toglabel(ci+1)=mode(nearborlabel(nearborlabel~=dlilabel));
    else
        toglabel(ci+1)=dlilabel;
    end
end
dli=find(gradientmean<GradTresh);
% find backroundground label
bl=mode(label(bw==255));
for di=1:length(dli)
    
    label(label==dli(di))=toglabel(dli(di));
end
RGBLabel=label2rgb(label);
figure,imshow(RGBLabel);title('segment result after valiation');
figure,plot(1:N,I(mingr,:),1:N,P_new(mingr,:));
legend('gray level','Treshold surface');
end