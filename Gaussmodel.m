function p= Gaussmodel(X,M,SIGMA)
d=-1/2;
 p = (((((2*pi)^(d))*(SIGMA))^(-(1/2)))   .*  exp((-1/2).*bsxfun(@minus,X,M).^2 ./(SIGMA)));
end