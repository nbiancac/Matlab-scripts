function Zleff=Zl_eff(freq,Z,machine,l,type,tol)
% Zl_eff(freq,Z,f0,Qs,l,'type')
clight = constants('clight');
Qs=machine.Qs;
taub=machine.taub;
f0=machine.f0;

            
sigma_z=taub*clight;

f_spec=[-1/taub:1/100/taub:1/taub]';
if freq(1)~=0
    f1=[-freq(end:-1:1);freq(1:end)];
    z1=[conj(Z(end:-1:1));Z(1:end)]; % NB only conjugation, no change in sign
elseif freq(1)==0
    f1=[-freq(end:-1:2);0;freq(2:end)];
    z1=[conj(Z(end:-1:2));Z(1:end),Z(2:end)]; % NB only conjugation, no change in sign
end
f=unique(sort([f1;f_spec]));
z=interp1(f1,z1,f);

omega=2*pi*f;
omega_0=2*pi*f0;
omega_s=Qs*omega_0;

if strcmp(type,'Gaussian')
    h=(omega*sigma_z/clight).^(2*l).*exp(-omega.^2*sigma_z^2/clight^2);
    
    nperc=tol;
    omega_part=omega(omega>0);
    h_part=h(omega>0);
    [max_h,max_ind]=max(h_part);
    omega_max=omega_part(max_ind);
    h_part2=h_part(omega_part>omega_max);
    omega_part2=omega_part(omega_part>omega_max);
    [h_part2,ind]=unique(h_part2); % delete trailed zeros;
    omega_part2=omega_part2(ind);
    omega_extr=interp1(h_part2,omega_part2,nperc/100*max_h);
    
    Pmax=floor(omega_extr/omega_0); % approx value for max intergation;
    p=-Pmax:Pmax;
    omega_p=p*omega_0+l*omega_s;

    omega_p=unique(omega_p);
    h_p=interp1(omega,h,omega_p);
    Z_p=interp1(omega,z,omega_p);
    N=(Z_p)./omega_p*omega_0.*h_p;
    D=h_p;
    
%     figure(1);
%     plot(f,h); xlim([-2e9 2e9]);
%     hold on; plot(f,imag(z)./omega*omega_0,'-r');
    
end
sum(D);
sum(N);
Zleff=(sum(N)./sum(D));
% figure(1)
% plot(omega_p,h_p*max(imag(Z_p)./omega_p)/max(h_p)); hold on;
% xlim([0 1e10])
% pause
end
    