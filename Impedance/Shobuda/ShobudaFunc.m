function [freq_sho,Zlong_sho,cond_vec]=ShobudaFunc(beta,a,t,L,material,fin,fstep,fout,fdiscrete,SaveDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               Shobuda.m
%   
%   Implementation of Shobuda's formulas for short thin insert (g<<a)
%   t=1e-3;        % thickness [m]
%   L=0.08;        % gap [m]
%   a=0.05;         % pipe radius [m]
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Longitudinal Impedance
Y1=[];Y2=[];Y3=[];
freq_sho=[];
Zlong_sho=[];

gamma=1/sqrt(1-beta^2);     % lorentz factor


Z0=120*pi;      % vacuum impedence
SS=200;           % saturation index for Ypole
js=besselzero(0,SS,1);
beta=sqrt(1-1/gamma^2);
sl=299792456.2;
f_vec=[];
cond_vec=[];

if ~isempty(fdiscrete)
    f_vec=fdiscrete;
    if size(fdiscrete,1)>size(fdiscrete,2) f_vec=f_vec'; end
else
    f_vec=fin:fstep:fout; 
end

for f=f_vec % frequency [Hz]
   
    omega=2*pi*f;
    K=omega/(beta*sl);
    Kr=K/gamma;
    bs=conj(sqrt(K^2*beta^2*a^2-js.^2));
    
    F     = 0*(K^2*beta^2*L/2)/(pi*sqrt(1i*K*beta*Z0*(material.sigma+1i*K*beta*eval(material.e_r)/Z0)))*(tanh( sqrt(1i*K*beta*Z0*(material.sigma+1i*K*beta*eval(material.e_r)/Z0))*t));
    Ypole = (- sum(     (4*pi*a)./(L*bs.^2) .* (1-exp(-1i*L.*bs/2/a))    ));
    Ycut  = (2*sqrt(2)*(1-1i)) / (sqrt(K*beta*L));
    G     = 0*(K^2*beta^2*L/2*Ypole*Ycut)/(pi*sqrt(1i*K*beta*Z0*(material.sigma+1i*K*beta*eval(material.e_r)/Z0)));
    Yres  = - ((  2*pi*sqrt(1i*K*beta*Z0*(material.sigma+1i*K*beta*eval(material.e_r)/Z0))  )/(K^2*beta^2*L) + G)*(tanh( sqrt(1i*K*beta*Z0*(material.sigma+1i*K*beta*eval(material.e_r)/Z0))*t));
    
    Z     = Z0*(1+F)/(1i*beta*K*a*besseli(0,Kr*a)^2)/(Ypole+Ycut+Yres);
    Zlong_sho = [Zlong_sho, Z];
    freq_sho  = [freq_sho, f];
    Y1=[Y1,Ypole];
    Y2=[Y2,Ycut];
    Y3=[Y3,Yres];
    cond_val=K^2*beta^2*(L/2)^2/sqrt(1i*K*beta*Z0*(material.sigma+1i*K*beta*eval(material.e_r)/Z0)*L/2)*...
         tanh(1i*K*beta*Z0*(material.sigma+1i*K*beta*eval(material.e_r)/Z0)*t);
    cond_vec=[cond_vec,cond_val];
end

%%
figure(1);
    subplot(211)
    plot(freq_sho,real(Zlong_sho),'k'); hold on;
    grid on;
    xlabel('f [Hz]')
    ylabel('Z_l [\Omega]')
    legend('Re(Z_l)')
    
    
    subplot(212)
    plot(freq_sho,imag(Zlong_sho),'b'); hold on;
    grid on;
    xlabel('f [Hz]')
    ylabel('Z_l [\Omega]')
    legend('Im(Z_l)');
%% Save
nome_fileRe=['Shobudalong_Re_G',num2str(L),'_Beta',num2str(beta),'_a',num2str(a),'_t',num2str(t),'_Material_',material.name,...
         '_fmin',num2str(min(f_vec)),...
         '_fstep',num2str(fstep),'_fmax',num2str(max(f_vec)),'.txt'];
nome_fileIm=['Shobudalong_Im_G',num2str(L),'_Beta',num2str(beta),'_a',num2str(a),'_t',num2str(t),'_Material_',material.name,...
         '_fmin',num2str(min(f_vec)),...
         '_fstep',num2str(fstep),'_fmax',num2str(max(f_vec)),'.txt'];  
% Write impedance
dlmwrite([SaveDir,nome_fileRe],[freq_sho',real(Zlong_sho')]);
dlmwrite([SaveDir,nome_fileIm],[freq_sho',imag(Zlong_sho')]);
disp(['Impedance data saved in: ',SaveDir]);
disp([nome_fileRe]);
disp([nome_fileIm]);
end