%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               Shobuda.m
%   
%   Implementation of Shobuda's formulas for short thin insert (g<<a)
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

set(0,'defaultlinelinewidth',2)
set(0,'defaultaxesfontname','timesnewroman');
%%
mainDir='/afs/cern.ch/user/n/nbiancac/scratch0/FiniteLength/Shobuda';
outputDir='/afs/cern.ch/user/n/nbiancac/scratch0/FiniteLength/Shobuda/Output';
cd(mainDir);
mkdir(outputDir);

%% Longitudinal Impedance
Y1=[];Y2=[];Y3=[];
freq_sho=[];
Zlong_sho=[];
t=1e-3;        % thickness [m]
g=0.08;        % gap [m]
a=0.05;         % pipe radius [m]
gamma=1000;     % lorentz factor
sigma_c=1e6; % conductivity [S/m]

Z0=120*pi;      % vacuum impedence
SS=20;           % saturation index for Ypole
js=besselzero(0,SS,1);
beta=sqrt(1-1/gamma^2);
sl=299792456.2;
eps_r=10;      % dielectric constant [F/m]

for f=1e8:1e6:1e10 % frequency [Hz]
   
    omega=2*pi*f;
    K=omega/(beta*sl);
    Kr=K/gamma;
    bs=conj(sqrt(K^2*beta^2*a^2-js.^2));
    
    F     = (K^2*beta^2*g/2)/(pi*sqrt(1i*K*beta*Z0*(sigma_c+1i*K*beta*eps_r/Z0)))*(tanh( sqrt(1i*K*beta*Z0*(sigma_c+1i*K*beta*eps_r/Z0))*t));
    Ypole = (- sum(     (4*pi*a)./(g*bs.^2) .* (1-exp(-1i*g.*bs/2/a))    ));
    Ycut  = (2*sqrt(2)*(1-1i)) / (sqrt(K*beta*g));
    G     = (K^2*beta^2*g/2*Ypole*Ycut)/(pi*sqrt(1i*K*beta*Z0*(sigma_c+1i*K*beta*eps_r/Z0)));
    Yres  = - ((  2*pi*sqrt(1i*K*beta*Z0*(sigma_c+1i*K*beta*eps_r/Z0))  )/(K^2*beta^2*g) + G)*(tanh( sqrt(1i*K*beta*Z0*(sigma_c+1i*K*beta*eps_r/Z0))*t));
    
    Z     = Z0*(1+F)/(1i*beta*K*a*besseli(0,Kr*a)^2)/(Ypole+Ycut+Yres);
    Zlong_sho = [Zlong_sho, Z];
    freq_sho  = [freq_sho, f];
    Y1=[Y1,Ypole];
    Y2=[Y2,Ycut];
    Y3=[Y3,Yres];
    
end

%%
figure(1);
    subplot(211)
    loglog(freq_sho,real(Zlong_sho),'k'); hold on
    grid on;
    xlabel('Frequency [Hz]')
    ylabel('Zlong [Ohm]')
    legend('Real(Z_l_o_n_g)')
    
    
    subplot(212)
    semilogx(freq_sho,imag(Zlong_sho),'b'); hold on;
    grid on;
    xlabel('Frequency [Hz]')
    ylabel('Zlong [Ohm]')
    legend('Imag(Z_l_o_n_g)');