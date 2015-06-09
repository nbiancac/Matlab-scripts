function Z3=junct(Z1,Z2,Freq);
% Calculate combined impedances 
%
% Usage : Z3=junct(Z1,Z2,Freq)

Z3=(Z1.*Z2)./(Z1+Z2);

