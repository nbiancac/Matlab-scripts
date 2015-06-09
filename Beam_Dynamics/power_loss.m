function [ploss,kloss]=power_loss(freq,Zlong,implength,machine,Pmax)
              
%%%%%%%%
% Calculates power loss [W] and loss factor in [Ohm/C]
% [ploss,kloss]=power_loss(freq,Zlong,implength,machine,Pmax)
% for a Gaussian spectrum
gamma=machine.gamma;
circum=machine.circ;
taub=machine.taub;
Nb=machine.Nb;
M=machine.M;
e=1.6021e-19;
mp=1.672623e-27;
clight=299792458;
beta=sqrt(1-1/gamma^2);
frev=beta*clight/circum;
sigz=taub*clight;

p=1:Pmax;
omega=2*pi*freq;
omega_rev=2*pi*frev;
omega_sample=omega_rev*M*p;
freq_sample=omega_sample/2/pi;
powspec=exp(-(omega.^2*taub^2));

Zlong_sample=interp1(omega,real(Zlong),omega_sample);
powspec_sample=interp1(omega,powspec,omega_sample);
% figure(3);
% loglog(freq,real(Zlong)); hold on;
% loglog(freq,powspec,'-r');
% loglog(freq_sample,real(Zlong_sample),'-ob'); hold on;
% loglog(freq_sample,powspec_sample,'-or'); hold off;
sum_impspec=sum(Zlong_sample.*powspec_sample);
% sum_impspec=trapz(freq,real(Zlong).*powspec)/M/frev;
ploss=2*(e*M*Nb*frev)^2*sum_impspec;
kloss=1/pi*trapz(omega,real(Zlong).*powspec);

end

