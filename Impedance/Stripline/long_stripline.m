function Z=long_stripline(freq,Zs,phi0,length,Ns)
% function Z=long_stripline(freq,Zs,phi0,length,Ns);

addpath(genpath('/afs/cern.ch/user/n/nbiancac/scratch0/Matlab-scripts/'));
c=constants('clight');
l=length;
omega=2*pi*freq;
Z=Ns*Zs*(phi0/2/pi)^2*(2*sin(omega*l/c).^2+1i*sin(2*omega*l/c));

end