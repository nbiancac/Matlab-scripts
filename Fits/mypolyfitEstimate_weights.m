function [estimate,uncertainty]=mypolyfitEstimate_weights(x,y,weigth,N)
% It gives the prediction and its uncertainty according to a polynomial fit 
% of order N of y data with no uncertainty
%
% x             x data
% y             y data
% N             order of the polynomial fit
% estimate      least square estimation of the y data points
% uncertainty   uncertainty of the least square estimation
%
% function [estimate,uncertainty]=mypolyfitEstimate(x,y,N)
% v. 1.0 Andrea Mostacci, 01-10-08

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

% weigth mus be a raw vector
[raw2,col2]=size(weigth);
if raw > col,
    weigth=weigth';
end

% Vectors inizialization
coeff=zeros(1,N+1);
unc_coeff=zeros(1,N+1);

Nsample=length(x);
% Weight matrix according to the book
WW= diag(1./(weigth.^2));
% Matrix computation according to the book
AA=zeros(length(x),N+1)+1;
for jjcol=1:N+1,
    AA(:,jjcol)=x.^(jjcol-1);
end
inverseMat=inv(AA.'*WW*AA);
%Equation 4.13
coeff=inverseMat*(AA.'*WW*y');
% The estimate as in eq. 2.26
estimate=AA*coeff;
% The uncertainty as in eq 4.17
uncertainty=sqrt(diag((AA*(inverseMat)'*AA.')));
