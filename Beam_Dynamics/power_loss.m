function [ploss,kloss]=power_loss(freq,Zlong,machine)
              
%%%%%%%%
% Calculates power loss [W] and loss factor in [Ohm/C]
% [ploss,kloss]=power_loss(freq,Zlong,machine)
% for a Gaussian spectrum up to 1% of magnitude of the spectrum  

gamma=machine.gamma;
circum=machine.circ;
taub=machine.taub; % 1 sigmaz in [s]
Nb=machine.Nb;
M=machine.M;
e=1.6021e-19;
mp=1.672623e-27;
clight=299792458;
beta=sqrt(1-1/gamma^2);
frev=beta*clight/circum;
sigz=taub*clight;

nperc=1;
omega=2*pi*freq;
omega_rev=2*pi*frev;
h=exp(-(omega.^2*taub^2));
    omega_part=omega(omega>0);
    h_part=h(omega>0);
    [max_h,max_ind]=max(h_part);
    omega_max=omega_part(max_ind);
    h_part2=h_part(omega_part>omega_max);
    omega_part2=omega_part(omega_part>omega_max);
    [h_part2,ind]=unique(h_part2); % delete trailed zeros;
    omega_part2=omega_part2(ind);
    omega_extr=interp1(h_part2,omega_part2,nperc/100*max_h);
    
    Pmax=floor(omega_extr/(omega_rev*M)); % approx value for max intergation;

p=1:Pmax;

omega_sample=omega_rev*M*p;
freq_sample=omega_sample/2/pi;
powspec=h;

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

