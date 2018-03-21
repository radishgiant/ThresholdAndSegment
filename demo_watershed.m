I=imread('cameraman.tif');
varn=EstimaOfNoise(I);
% [SI,SNR]=ImageReduNoise(I,varn,'NeigborHoodSize',[7,7],'SigLevel','SigLevel3');
Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);
figure, imshow(Iobr), title('Opening-by-reconstruction (Iobr)')
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
SI = imcomplement(Iobrcbr);
figure, imshow(SI), title('Opening-closing by reconstruction (Iobrcbr)')
% gradient Image
h1=[-1,0,1];
h2=h1';
Ix=imfilter(double(SI),h1,'replicate');
Iy=imfilter(double(SI),h2,'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);%ÇóÃþ
figure, imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
% candidate egde pixels that have the
% regional maxima 
CanI=imregionalmax(gradmag);
% smooth gradient Image Gs
h=fspecial('average',[3 3]);
Gs=imfilter(gradmag,h);
Gs(CanI~=0)=gradmag(CanI~=0);
bw=zeros(size(I));
T=4*sqrt(varn);
bw(Gs>T)=gradmag(Gs>T);
figure,imshow(bw,[]);title('bw')
label=watershed(bw);
RGBLabel=label2rgb(label);
figure,imshow(RGBLabel);
title(' watershed with treshold Image of gradient Image');
Rlabel=RegionMerging(I,label,1000,2);
RGBRlabel=label2rgb(Rlabel);
figure,imshow(RGBRlabel);title('watershed segment after RegionMerging');