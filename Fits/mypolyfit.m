function [coeff,unc_coeff]=mypolyfit(x,y,N)
% Polinomial fit of order N of data with no uncertainty on the y data
%
% x           x data
% y           y data 
% N           order of the polynomial fit
% coeff       coefficient of the least square estimation 
% unc_coeff   uncertainty of the coefficients
%
% y=coeff(1)*x^N+ coeff(2)*x^(N-1)+ ... + coeff(N);
%
% function [coeff,unc_coeff]=mypolyfit(x,y,N)
% v. 1.1 Andrea Mostacci, 01-10-08

% The routine gives the same coefficents of the polyfit routine but it
% seems faster. Uncertainty on the polynomial fit parameter is inferred 
% by the difference between the y value and its best estimate and it is 
% identical to regression2.m uncertainty. 
% Reference:    Least Squares Methods by V. Blobel
%               in Formulae and Methods in experimenta data evaluation
%                  with emphasis on High Energy Physics


if length(x)~= length(y)
    error('Vectors must be of the same length');
end

% y mus be a raw vector
[raw,col]=size(y);
if raw > col,
    y=y';
end

% Vectors inizialization
coeff=zeros(1,N+1);
unc_coeff=zeros(1,N+1);

Nsample=length(x);
% Matrix computation according to the book
AA=zeros(length(x),N+1)+1;
for jjcol=1:N+1,
    AA(:,jjcol)=x.^(jjcol-1);
end
inverseMat=inv(AA.'*AA);
coeff=inverseMat*(AA.'*y');
% To estimate the variance, we compute the sum of square residuals
% as in eq. 2.24
SquaredSumResiduals=(y*y'-coeff'*AA.'*y')/(Nsample-N);
% The uncertainty as in eq 2.21
unc_coeff=sqrt(SquaredSumResiduals*diag(inverseMat));
% To be coherent with the polyfit syntax, we invert the order of the vector
index=N+1:-1:1;
coeff=coeff(index);
unc_coeff=unc_coeff(index);
