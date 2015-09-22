function [Zeff,kick_f,omega_shifted,mat]=Zt_eff(freq,Z,machine,m,nx,type,particle,plane,spec,nperc)
% function [Zeff,kick_f,omega_shifted,mat]=Zt_eff(freq,Z,machine,m,nx,type,particle,plane,spec,nperc)
% requires changing machine.taub (not sensitive on machine.sigmaz that is recalculated)

[e,m0,clight,E0]=particle_param(particle);
taub=machine.taub;
tunes=machine.Qs;
tune=eval(['machine.Q',plane]);
chroma=eval(['machine.chroma',plane]);
gammatr=machine.gammatr;
gammarel=machine.gamma;
f0=machine.f0;
T0=1/f0;
Qs=machine.Qs;
beta=machine.beta;
Nb=machine.Nb;
Lb=4*taub*beta*clight;
Ib=Nb*e*f0;
sigma_z=taub*beta*clight;
M=machine.M;
r0=e^2/(m0*clight^2);

if size(Z,1)<size(Z,2)
    Z=conj(Z');
end
if size(freq,1)<size(freq,2)
    freq=(freq)';
end
f_spec=[-1/taub:1/100/taub:1/taub]';
if freq(1)~=0
    f1=[-freq(end:-1:1);freq(1:end)];
    z1=[-conj(Z(end:-1:1));Z(1:end)]; % NB conjugation and change in sign
elseif freq(1)==0
    f1=[-freq(end:-1:2);0;freq(2:end)];
    z1=[-conj(Z(end:-1:2));Z(1);Z(2:end)]; % NB conjugation and change in sign
end
f=unique(sort([f1;f_spec]));
omega=unique(2*pi*f);
f=omega/2/pi;
z=interp1(f1,z1,f,'linear','extrap');


omega_0=2*pi*f0;
omega_s=Qs*omega_0;
eta=1/gammatr^2-1/gammarel^2;
omega_beta=omega_0*tune;
omega_xi=omega_beta*chroma/eta;


if strcmp(type,'Gaussian')
    h=((omega-omega_xi)*sigma_z/clight).^(2*m).*exp(-(omega-omega_xi).^2*sigma_z^2/clight^2);
    
    [ind0,val]=find(h==0);
    [indnan,val]=find(isnan(h)==1);
    
    inddel=[ind0;indnan];
    omega(inddel)=[];    f(inddel)=[];    z(inddel)=[];    h(inddel)=[];
    
    if strcmp(spec,'h_range')
        omega_part=omega(omega>omega_xi);
        h_part=h(omega>omega_xi);
        [max_h,max_ind]=max(h_part);
        omega_max=omega_part(max_ind);
        h_part2=h_part(omega_part>omega_max);
        omega_part2=omega_part(omega_part>omega_max);
        [h_part2,ind]=unique(h_part2); % delete trailed zeros;
        omega_part2=omega_part2(ind);
        omega_extr=interp1(h_part2,omega_part2,nperc/100*max_h);
        f_extr=omega_extr/2/pi;
        Pmax=floor((omega_extr-(nx+tune+m*tunes)*omega_0)/(M*omega_0)); % approx value for max intergation;
        disp(['Pmax=',num2str(Pmax)])
        if Pmax>1e6; error('max Pmax is > 1e6'); end
        p_prime=-Pmax:Pmax;
        omega_p=(nx+p_prime*M+tune)*omega_0+m*omega_s;   
    elseif strcmp(spec,'z_range')
        omega_part=omega(omega>0);
        z_part=real(z(omega>0));
        [max_z,max_ind]=max(z_part);
        omega_max=omega_part(max_ind);
        
        z_part2=z_part(omega_part>=omega_max);
        omega_part2=omega_part(omega_part>=omega_max);
        [z_part2,ind]=unique(z_part2); % delete trailed zeros;
        omega_part2=omega_part2(ind);
        omega_2=interp1(z_part2,omega_part2,nperc/100*max_z);
        
        z_part1=z_part(omega_part<=omega_max);
        omega_part1=omega_part(omega_part<=omega_max);
        [z_part1,ind]=unique(z_part1); % delete trailed zeros;
        omega_part1=omega_part1(ind);
        omega_1=interp1(z_part1,omega_part1,nperc/100*max_z);
        
        Pmax=floor((omega_2-(nx+tune+m*tunes)*omega_0)/(M*omega_0)); % approx value for max intergation;
        Pmin=floor((omega_1-(nx+tune+m*tunes)*omega_0)/(M*omega_0)); % approx value for max intergation;
%         disp(['Pmax=',num2str(abs(Pmax)),', Pmin=',num2str(abs(Pmin))])
        if Pmax>1e6; error('Pmax is > 1e6'); end
        p_prime=([-Pmax:-Pmin,Pmin:Pmax]);
        omega_p=(nx+p_prime*M+tune)*omega_0+m*omega_s;
    end
    
    kick_f=trapz(f,z.*h);    
    h_p=interp1(omega,h,omega_p);
    Z_p=interp1(omega,z,omega_p);
    N=(Z_p).*h_p;
    
    %   

    h=((omega-omega_xi)*sigma_z/clight).^(2*m).*exp(-(omega-omega_xi).^2*sigma_z^2/clight^2);
    
    [ind0,val]=find(h==0);
    [indnan,val]=find(isnan(h)==1);
    inddel=[ind0;indnan];
    omega(inddel)=[];    f(inddel)=[];    z(inddel)=[];    h(inddel)=[];
    
    omega_part=omega(omega>omega_xi);
    h_part=h(omega>omega_xi);
    [max_h,max_ind]=max(h_part);
    omega_max=omega_part(max_ind);
    h_part2=h_part(omega_part>omega_max);
    omega_part2=omega_part(omega_part>omega_max);
    [h_part2,ind]=unique(h_part2);
    omega_part2=omega_part2(ind);
    omega_extr=interp1(h_part2,omega_part2,nperc/100*max_h);
    f_extr=omega_extr/2/pi;
    Pmax=floor((omega_extr-(nx+tune+m*tunes)*omega_0)/(M*omega_0)); % approx value for max intergation;
    disp(['Pmax=',num2str(Pmax)])
    if Pmax>1e6; error('max Pmax is > 1e6'); end
    p_prime=-Pmax:Pmax;
    omega_p2=(nx+p_prime*M+tune)*omega_0+m*omega_s;
    h_p2=interp1(omega,h,omega_p2);
    D=h_p2;
%     D=h_p;
    
    kick_f=trapz(f,z.*h);
    Zeff=(sum(N)./sum(D));
    omega_shifted=1/4/pi*gamma(m+1/2)/(2^m*factorial(m))*(machine.Nb*r0*clight^2)/(gammarel*T0*tune*omega_0*sigma_z)*1i*Zeff;

elseif strcmp(type,'Sinusoidal')
    
    h=(4*taub)^2/(2*pi^4)*(abs(m)+1)^2*(1+(-1)^abs(m).*cos((omega-omega_xi)*(4*taub)))./(((omega-omega_xi)*(4*taub)/pi).^2-(abs(m)+1)^2).^2;
    [ind0,val]=find(h==0);
    [indnan,val]=find(isnan(h)==1);
    
    inddel=[ind0;indnan];
    omega(inddel)=[];    f(inddel)=[];    z(inddel)=[];    h(inddel)=[];
    
    if strcmp(spec,'h_range')
        omega_part=omega(omega>omega_xi);
        h_part=h(omega>omega_xi);
        [max_h,max_ind]=max(h_part);
        omega_max=omega_part(max_ind);
        h_part2=h_part(omega_part>omega_max);
        omega_part2=omega_part(omega_part>omega_max);
        [h_part2,ind]=unique(h_part2); % delete trailed zeros;
        omega_part2=omega_part2(ind);
        omega_extr=interp1(h_part2,omega_part2,nperc/100*max_h);
        f_extr=omega_extr/2/pi;
        Pmax=floor((omega_extr-(nx+tune+m*tunes)*omega_0)/(M*omega_0)); % approx value for max intergation;
%         disp(['Pmax=',num2str(Pmax)])
        if Pmax>1e6; error('max Pmax is > 1e6'); end
        p_prime=-Pmax:Pmax;
        omega_p=(nx+p_prime*M+tune)*omega_0+m*omega_s;   
    elseif strcmp(spec,'z_range')
        omega_part=omega(omega>0);
        z_part=real(z(omega>0));
        [max_z,max_ind]=max(z_part);
        omega_max=omega_part(max_ind);
        
        z_part2=z_part(omega_part>=omega_max);
        omega_part2=omega_part(omega_part>=omega_max);
        [z_part2,ind]=unique(z_part2); % delete trailed zeros;
        omega_part2=omega_part2(ind);
        omega_2=interp1(z_part2,omega_part2,nperc/100*max_z);
        
        z_part1=z_part(omega_part<=omega_max);
        omega_part1=omega_part(omega_part<=omega_max);
        [z_part1,ind]=unique(z_part1); % delete trailed zeros;
        omega_part1=omega_part1(ind);
        omega_1=interp1(z_part1,omega_part1,nperc/100*max_z);
        
        Pmax=floor((omega_2-(nx+tune+m*tunes)*omega_0)/(M*omega_0)); % approx value for max intergation;
        Pmin=floor((omega_1-(nx+tune+m*tunes)*omega_0)/(M*omega_0)); % approx value for max intergation;
        disp(['Pmax=',num2str(abs(Pmax)),', Pmin=',num2str(abs(Pmin))])
        disp(['f_max=',num2str((omega_2/2/pi)),', f_min=',num2str((omega_1/2/pi))])
        if Pmax>1e6; error('Pmax is > 1e6'); end
        p_prime=([-Pmax:-Pmin,Pmin:Pmax]);
        omega_p=(nx+p_prime*M+tune)*omega_0+m*omega_s;
    end
    
    kick_f=trapz(f,z.*h);    
    h_p=interp1(omega,h,omega_p);
    Z_p=interp1(omega,z,omega_p);
    N=(Z_p).*h_p;
%   

    h=(4*taub)^2/(2*pi^4)*(abs(m)+1)^2*(1+(-1)^abs(m).*cos((omega-omega_xi)*(4*taub)))./(((omega-omega_xi)*(4*taub)/pi).^2-(abs(m)+1)^2).^2;
    [ind0,val]=find(h==0);
    [indnan,val]=find(isnan(h)==1);
    inddel=[ind0;indnan];
    omega(inddel)=[];    f(inddel)=[];    z(inddel)=[];    h(inddel)=[];
    
    omega_part=omega(omega>omega_xi);
    h_part=h(omega>omega_xi);
    [max_h,max_ind]=max(h_part);
    omega_max=omega_part(max_ind);
    h_part2=h_part(omega_part>omega_max);
    omega_part2=omega_part(omega_part>omega_max);
    [h_part2,ind]=unique(h_part2);
    omega_part2=omega_part2(ind);
    omega_extr=interp1(h_part2,omega_part2,nperc/100*max_h);
    f_extr=omega_extr/2/pi;
    Pmax=floor((omega_extr-(nx+tune+m*tunes)*omega_0)/(M*omega_0)); % approx value for max intergation;
    disp(['Pmax=',num2str(Pmax)])
    if Pmax>1e6; error('max Pmax is > 1e6'); end
    p_prime=-Pmax:Pmax;
    omega_p2=(nx+p_prime*M+tune)*omega_0+m*omega_s;
    h_p2=interp1(omega,h,omega_p2);
    D=h_p2;
%     D=h_p;
    Zeff=(sum(N)./sum(D));
    omega_shifted=(abs(m)+1)^-1*(1i*e*beta*Ib)/(2*m0*gammarel*tune*omega_0*Lb)*Zeff; 
end

mat=[omega_p'/2/pi,conj(Z_p'),h_p'];
% mat=[f,z,h];

end
    