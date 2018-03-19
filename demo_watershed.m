I=imread('cameraman.tif');
varn=EstimaOfNoise(I);
SI=ImageReduNoise(I,varn,'NeigborHoodSize',[13,13],'SigLevel','SigLevel3');
figure,imshow(SI);title('smoothed Image');
T=sqrt(varn);

% gradient Image
h1=[-1,0,1];
h2=h1';
Ix=imfilter(double(SI),h1,'replicate');
Iy=imfilter(double(SI),h2,'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);%ÇóÃþ
figure,imshow(gradmag./max(gradmag(:)));title('gradient Image');
% candidate egde pixels
h=fspecial('laplacian',0);
CanI=imfilter(I,h);
% smooth gradient Image Gs
h=fspecial('average',[3 3]);
Gs=imfilter(gradmag,h);
Gs(CanI~=0)=gradmag(CanI~=0);
bw=zeros(size(I));

T=10;
bw(Gs>T)=gradmag(Gs>T);
figure,imshow(bw);title('bw')
label=watershed(bw);
RGBLabel=label2rgb(label);
figure,imshow(RGBLabel);
title(' watershed with treshold Image of gradient Image');
Rlabel=RegionMerging(I,label,1000,5);
RGBRlabel=label2rgb(Rlabel);
figure,imshow(RGBRlabel);title('watershed segment after RegionMerging');