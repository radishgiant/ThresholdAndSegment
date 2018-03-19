function varn=EstimaOfNoise(NI)
% 图像转为double型和噪声叠加
h1=fspecial('laplacian',0);%滤波器L1
h2=fspecial('laplacian',1);%滤波器L2
N=(h2-h1)*2;
imI=imfilter(NI,N);
%估计噪声方差
varn=sum(abs(imI(:)))./(size(imI,1)-2)./(size(imI,2)-2)./6.*sqrt(pi/2);
varn=varn.^2;
end