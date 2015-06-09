function kloss=kl_factor(freq,Z,taub, f0,Qs,l,type)
% Zl_eff(freq,Z,f0,Qs,l,'type')
clight = 299792458;

            
sigma_z=taub*clight;

f_spec=[-1/taub:1/100/taub:1/taub]';
f1=[-freq(end:-1:1);freq(2:end)];
z1=[Z(end:-1:1);Z(2:end)];
f=unique(sort([f1;f_spec]));
z=interp1(f1,z1,f);

omega=2*pi*f;
omega_0=2*pi*f0;
omega_s=Qs*omega_0;

if strcmp(type,'Gaussian')
    h=(omega*sigma_z/clight).^(2*l).*exp(-omega.^2*sigma_z^2/clight^2);
    
    kloss=trapz(f,z.*h);
    
end
% figure(1)
% plot(omega_p,h_p*max(imag(Z_p)./omega_p)/max(h_p)); hold on;
% xlim([0 1e10])
% pause
end
    