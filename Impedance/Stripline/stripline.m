function [Zl,Zdip]=stripline(freq,Zs,phi0,length,R,Ns)
% function Z=long_stripline(freq,Zs,phi0,length,R,Ns);

addpath(genpath('/afs/cern.ch/user/n/nbiancac/scratch0/Matlab-scripts/'));
c=constants('clight');
l=length;
omega=2*pi*freq;
Zl=Ns*Zs*(phi0/2/pi)^2*(2*sin(omega*l/c).^2+1i*sin(2*omega*l/c));
Zdip=Zl./(2*pi*freq)*c/R^2*(4/phi0)^2*sin(phi0/2)^2;
end