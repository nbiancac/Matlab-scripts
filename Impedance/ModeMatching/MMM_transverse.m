%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               MMM_transverse.m
%   
%   Impedance of a Finite length device with the Mode Matching Method (v2.05.2014). 
%
%   function [freq, Z, TMmodes_WG, TEmodes_WG, TMmodes_C, TEmodes_C] = MMM_transverse(beta,b,t,L,material,fin,fstep,fout,fdiscrete,f_post, P,S,SaveDir,fields, ph, loss)
%   
%   Function description:
%     This function calculates the dipolar (or driving) coupling impedance for a
%     cylindrical insert of length L, inner radius b and thickness t, surrounded by PEC
%     and accessed with two beam pipes of radius b. The insert is made of a
%     linear isotropic homogeneus dispersive material specified by its
%     permeability mu(omega) and permittivity eps(omega) complex quantities. 
%   
%   Input:
%     beta      :       relativistic beta,
%     b         :       inner radius               [m]
%     t         :       thickness of the insert    [m]
%     L         :       length of the insert       [m]
%     material  :       it is a structure that contains 4 fields:
%           ".e_r       :        relative permettivity    
%           ".mu_r      :        relative permeability     
%           ".sigma     :        conductivity      [S/m]
%           NB: mu_r can be either an expression mu_r(f) to be evaluated
%           in function of the frequency "f" or a table for a given
%           frequency set to be specified in fdiscrete. Similarly holds for
%           e_r.
%     fin       :       initial frequency          [Hz]
%     fstep     :       step frequency             [Hz]
%     fout      :       final frequency            [Hz]
%     fdiscrete :       if given, the code will use this set of frequencies, otherwise you can 
%                       just put [] and it will be ignored. 
%     P         :       number of transverse modes
%     S         :       number of longitudinal modes
%     f_post    :       frequencies in which you want to postprocess. Until
%                       now only Power loss calculation in the material is done.
%     SaveDir   :       Directory where you save impedances and figures.
%     fields    :       Calculates fields pattern at each frequency in f_post on the plane cut "plane";
%     ph        :       Used in the fields calculation, indicates the angle
%                       of the plane with respect to the x plane. X plane is ph=0, Y plane is ph=pi/2.
%     loss      :       Not yet implemented! It will calculate power loss at each frequency in  f_post;
%
%   Output:
%     freq      :       Frequencies used [Hz]
%     Zdip     :        Dipolar impedance [Ohm] 
%     TM(E)modes_WG:       TM(E) modes frequencies in the pipes [GHz]  
%     TM(E)modes_C :       TM(E) modes frequencies in the empty cavity [GHz]
%  
%   Written Files:
%     "MMMdip_Re_*.txt": Real part of Zdip;
%     "MMMdip_Im_*.txt": Imaginary part of Zdip;
%     
% __________ Examples _________
% MMMdir='myWorkingDirectory';
% cd(MMMdir)
% beta=1;
% sl=299792448;
% b=0.05;
% t=0.04;
% L=0.0008;
% c=b+t;
% 
% material={};
% material.name='Ferrite4A4';
% mui=460;
% tauf=1/(2*pi*20e6);
% material.mu_r=['1+',num2str(mui),'/(1+1i*2*pi*f*',num2str(tauf),')'];
% material.e_r='12';
% material.sigma=1e-6;
% 
% 
% fin=1e8;
% fstep=1e6;
% fout=9e8;
% fdiscrete=[];
% f_post=[];
% 
% P=15;
% S=15;
% SaveDir='mySaveDirectory';
% [freq,Zdip,~,~]=MMM_transverse(beta,b,t,L,material,fin,fstep,fout,fdiscrete,f_post,P,S,SaveDir,0,0,0);
%     
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [f_vec, Z, TMmodes_WG, TEmodes_WG, TMmodes_C, TEmodes_C]=MMM_transverse(beta,b,t,L,material,fin,fstep,fout,fdiscrete,f_post, P,S,SaveDir,fields, ph,~)
    
% Beam

    Q=1; % Beam charge [C]
    gamma=1/sqrt(1-beta^2); % Relativistic gamma

% Physical constants

    sl=299792458; % speed of light [m/s]
    epsilon0=10^7/(4*pi*sl^2); % Vacuum dielectric constant [Farad/m]
    mu0=4*pi*10^-7; % Vacuum permeability [Henry^-1/m]
    Z0=sqrt(mu0/epsilon0);
    Y0=1/Z0;

% Geometry

    c=b+t;

% Material

    sigma=material.sigma;


% Frequency scan

    Zrew=[];
   
    if ~isempty(fdiscrete)
    f_vec=fdiscrete;
    if size(fdiscrete,1)>size(fdiscrete,2) fdiscrete=fdiscrete'; f_vec=fdiscrete; end
    else
        f_vec=fin:fstep:fout; 
    end
    if size(f_post,1)>size(f_post,2) f_post=f_post'; end
    f_vec=unique(sort([f_vec,fdiscrete,f_post]));
    
    nu=1;
%     if P>300; error('P_max=300 exceeded.'); end; 
%     [Jzero,Jzeroprime]=besselzeros2(nu+1,P);
%     Jzeroprime=Jzeroprime(:,2);
%     Jzero=Jzero(:,2);
    Jzero=zerobess('J',nu,P);
    Jzeroprime=zerobess('DJ',nu,P);
    
% Matrix inizialization

    B =(zeros(S,1));
    F1=(zeros(S,S));
    F2=(zeros(S,S));
    F3=(zeros(S,S));
    F4=(zeros(S,S));
    G1=(zeros(S,P));
    G2=(zeros(S,P));
    G3=(zeros(S,P));
    G4=(zeros(S,P));
    G5=(zeros(S,P));
    IIss=(zeros(S,S));
    M1=(zeros(P,S));
    M2=(zeros(P,S));
    M3=(zeros(P,S));
    M4=(zeros(P,S));
    M5=(zeros(P,S));
    N1=(zeros(P,P));
    N2=(zeros(P,P));
    N3=(zeros(P,P));
    T1=(zeros(S,S));
    T2=(zeros(S,S));
    T3=(zeros(S,S));
    T4=(zeros(S,S));
    T5=(zeros(S,S));
    FF1=(zeros(S,S));
    FF2=(zeros(S,S));
    UU1=(zeros(S,S));
    UU2=(zeros(S,S));

    Z_integral_TM_sx=(zeros(1,P));
    Z_integral_TM_dx=(zeros(1,P));
    Z_integral_TE_sx=(zeros(1,P));
    Z_integral_TE_dx=(zeros(1,P));

    Z_CTEirr=(zeros(1,P));
    Z_DTEirr=(zeros(1,P));
    Z_ATMirr=(zeros(1,S));
    Z_ATEirr=(zeros(1,S));

    Z_CTMsol=(zeros(1,P));
    Z_DTMsol=(zeros(1,P));
    Z_ATM1sol=(zeros(1,S));
    Z_CTEsol=(zeros(1,P));
    Z_DTEsol=(zeros(1,P));
    Z_ATM2sol=(zeros(1,S));
    Z_ATEsol=(zeros(1,S));
    Z=(zeros(1,length(f_vec)));

    A_TM_post=[];
    C_TM_post=[];
    D_TM_post=[];
    A_TE_post=[];
    C_TE_post=[];
    D_TE_post=[];  
    
for ff=1:length(f_vec)
                             
    f=f_vec(ff);
    disp(f);
    
    % Material properties
    
    omega=2*pi*f;
    if ischar(material.e_r)
        epsilonf=epsilon0*(eval(material.e_r)-1i*material.sigma/(omega*epsilon0)); % Relative permeability
    elseif length(material.e_r)==length(f_vec)
        epsilonf=epsilon0*((material.e_r(ff))-1i*material.sigma/(omega*epsilon0)); 
    end
    if ischar(material.mu_r)    % Relative susceptibility
        muf=mu0*eval(material.mu_r);
    elseif length(material.mu_r)==length(f_vec)
        muf=mu0*material.mu_r(ff);
    end
    Zf=sqrt(muf/epsilonf);
    Yf=1/Zf;
   

    % Kappas
    k0=omega/sl;
    kf=omega*sqrt(muf*epsilonf);
    kb=omega/(beta*sl);

    for s=0:S-1;
        
        ss=s+1; %indice per riempire le matrici
        if s==0; es=1;  else es=2;    end
        if s==0; delta_s0=1; else delta_s0=0; end
        
        for p=1:P;
            
            pp=p; % indici per riempire le matrici

            alpha_0           =       k0*b;
            alpha_02          =       alpha_0^2;
            alpha_f           =       kf*b;
            alpha_f2          =       alpha_f^2;
            alpha_b           =       kb*b;
            alpha_b2          =       alpha_b^2;
            alpha_p           =       Jzero(p);
            alpha_s           =       s*pi/(L/b);
            alpha_s2          =       alpha_s^2;
            alpha_p2          =       alpha_p^2;
            alpha_ps          =       sqrt(alpha_p^2+alpha_s^2);
            alpha_ps2         =       alpha_ps^2;
            alpha_tilde_s2    =       alpha_0^2-alpha_s^2;
            alpha_tilde_p2    =       alpha_0^2-alpha_p^2;
            alpha_tilde_sf2   =       alpha_f^2-alpha_s^2;
            alpha_tilde_sf    =       sqrt(alpha_tilde_sf2);
            alpha_tilde_p     =       conj(sqrt(alpha_tilde_p2));
            alpha_tilde_s     =       (sqrt(alpha_tilde_s2));
            
            beta_s            =       alpha_s;
            beta_s2           =       beta_s^2;
            beta_p            =       Jzeroprime(p);
            beta_p2           =       beta_p^2;
            beta_ps2          =       beta_p2+beta_s2;
            beta_tilde_p2     =       alpha_02-beta_p2;
            beta_tilde_s2     =       alpha_02-beta_s2;
            beta_tilde_p      =       conj(sqrt(beta_tilde_p2));
            beta_tilde_s      =       (sqrt(beta_tilde_s2));
            beta_ps           =       sqrt(beta_ps2);
           
            if sigma==0 
                alpha_tilde_sf    =       conj(alpha_tilde_sf);
            end
            
            in                 =        alpha_tilde_sf;
            W_TE_sf            =        ((besselh(0,1,(c*alpha_tilde_sf)/b,1) - besselh(2,1,(c*alpha_tilde_sf)/b,1))*besselh(1,2,alpha_tilde_sf,1) - exp(((2*1i)*(b - c)*alpha_tilde_sf)/b)*besselh(1,1,alpha_tilde_sf,1)*(besselh(0,2,(c*alpha_tilde_sf)/b,1) - besselh(2,2,(c*alpha_tilde_sf)/b,1)))/((besselh(0,1,(c*alpha_tilde_sf)/b,1) - besselh(2,1,(c*alpha_tilde_sf)/b,1))*besselh(1,2,alpha_tilde_sf,1));
            W_TM_sf            =        1 - (exp(((2*1i)*(b - c)*alpha_tilde_sf)/b)*besselh(1,1,alpha_tilde_sf,1)*besselh(1,2,(c*alpha_tilde_sf)/b,1))/(besselh(1,1,(c*alpha_tilde_sf)/b,1)*besselh(1,2,alpha_tilde_sf,1));
            W_TEprime_sf       =        (besselh(0,1,(c*alpha_tilde_sf)/b,1)*(besselh(0,2,alpha_tilde_sf,1) - besselh(2,2,alpha_tilde_sf,1)) + besselh(2,1,(c*alpha_tilde_sf)/b,1)*(-besselh(0,2,alpha_tilde_sf,1) + besselh(2,2,alpha_tilde_sf,1)) - exp(((2*1i)*(b - c)*alpha_tilde_sf)/b)*(besselh(0,1,alpha_tilde_sf,1) - besselh(2,1,alpha_tilde_sf,1))*(besselh(0,2,(c*alpha_tilde_sf)/b,1) - besselh(2,2,(c*alpha_tilde_sf)/b,1)))/(2*(besselh(0,1,(c*alpha_tilde_sf)/b,1) - besselh(2,1,(c*alpha_tilde_sf)/b,1))*besselh(1,2,alpha_tilde_sf,1));                  
            W_TMprime_sf       =        (-(exp((2*1i)*(b-c)/b*alpha_tilde_sf)*(besselh(0,1,alpha_tilde_sf,1) - besselh(2,1,alpha_tilde_sf,1))*besselh(1,2,(c*alpha_tilde_sf)/b,1)) + besselh(1,1,(c*alpha_tilde_sf)/b,1)*(besselh(0,2,alpha_tilde_sf,1) - besselh(2,2,alpha_tilde_sf,1)))/(2*besselh(1,1,(c*alpha_tilde_sf)/b,1)*besselh(1,2,alpha_tilde_sf,1));

            
            % Matrices for coefficients
            B(ss,1)            =     ((-1i/2)*alpha_b)/(b*pi*0.5*(besseli(0, alpha_b/gamma)-besseli(2, alpha_b/gamma))*(alpha_b^2 - alpha_s^2)) + ((1i/2)*(-1)^s*alpha_b)/(b*exp((1i*L*alpha_b)/b)*pi*0.5*(besseli(0, alpha_b/gamma)-besseli(2, alpha_b/gamma))*(alpha_b^2 - alpha_s^2));
            IIss(ss,ss)        =     (-1)^s;
            
            F1(ss,ss)          =     -((sqrt(L)*alpha_s*W_TE_sf)/(sqrt(2)*alpha_tilde_sf^2)) + (sqrt(L)*(besseli(0, alpha_s,1) - besseli(2, alpha_s,1))*alpha_f*alpha_s*W_TEprime_sf)/(sqrt(2)*(besseli(0, alpha_s,1) + besseli(2, alpha_s,1))*Yf*Z0*alpha_0*alpha_tilde_sf) + (sqrt(2)*sqrt(L)*(besseli(2, alpha_s,1)*besselj(0, alpha_tilde_s,1) + besseli(0, alpha_s,1)*besselj(2, alpha_tilde_s,1))*alpha_f*alpha_s*W_TEprime_sf)/((besseli(0, alpha_s,1) + besseli(2, alpha_s,1))*(besselj(0, alpha_tilde_s,1) - besselj(2, alpha_tilde_s,1))*Yf*Z0*alpha_0*alpha_tilde_sf) - (sqrt(2)*sqrt(L)*(besseli(2, alpha_s,1)*besselj(0, alpha_tilde_s,1) + besseli(0, alpha_s,1)*besselj(2, alpha_tilde_s,1))*alpha_f*alpha_s*delta_s0*W_TEprime_sf)/((besseli(0, alpha_s,1) + besseli(2, alpha_s,1))*(besselj(0, alpha_tilde_s,1) - besselj(2, alpha_tilde_s,1))*Yf*Z0*alpha_0*alpha_tilde_sf);
            F2(ss,ss)          =     ((-1i/2)*sqrt(L)*(besseli(0, alpha_s,1) - besseli(2, alpha_s,1))*sqrt(es)*W_TM_sf)/((besseli(0, alpha_s,1) + besseli(2, alpha_s,1))*Z0*alpha_0) - ((1i/2)*sqrt(L)*(besseli(0, alpha_s,1) - besseli(2, alpha_s,1))*delta_s0*sqrt(es)*W_TM_sf)/((besseli(0, alpha_s,1) + besseli(2, alpha_s,1))*Z0*alpha_0) + (1i*besselj(2, alpha_tilde_s,1)*alpha_0*sqrt(L/es)*W_TM_sf)/(besselj(1, alpha_tilde_s,1)*Z0*alpha_tilde_s) - ((1i/2)*sqrt(L)*(besseli(0, alpha_s,1) - besseli(2, alpha_s,1))*alpha_s^2*sqrt(es)*W_TM_sf)/((besseli(0, alpha_s,1) + besseli(2, alpha_s,1))*Z0*alpha_0*alpha_tilde_sf^2) - ((1i/2)*sqrt(L)*alpha_f^2*alpha_0*((alpha_0^2 - alpha_s^2)^(-1) + ((-besseli(0, alpha_s,1) + besseli(2, alpha_s,1))/(besseli(0, alpha_s,1) + besseli(2, alpha_s,1)) - ((besselj(0, alpha_tilde_s,1) + besselj(2, alpha_tilde_s,1))*alpha_s^2)/((besselj(0, alpha_tilde_s,1) - besselj(2, alpha_tilde_s,1))*(alpha_0^2 - alpha_s^2)))/alpha_0^2)*sqrt(es)*W_TM_sf)/(Z0*alpha_tilde_sf^2) + ((1i/2)*sqrt(L)*besselj(2, alpha_tilde_s,1)*alpha_0*alpha_s^2*sqrt(es)*W_TM_sf)/(Z0*alpha_tilde_s*(besselj(1, alpha_tilde_s,1) - besselj(0, alpha_tilde_s,1)*alpha_tilde_s)*alpha_tilde_sf^2) + (1i*alpha_f*sqrt(L/es)*W_TMprime_sf)/(Zf*alpha_tilde_sf);
            F3(ss,ss)          =     1 + ((besseli(0, alpha_s,1) - besseli(2, alpha_s,1))*alpha_f*alpha_s^2*sqrt(es)*W_TEprime_sf)/(sqrt(2)*(besseli(0, alpha_s,1) + besseli(2, alpha_s,1))*Yf*Z0*alpha_0*alpha_tilde_sf*W_TE_sf) - (alpha_f*(1-delta_s0)*(((besseli(0, alpha_s,1) - besseli(2, alpha_s,1))*alpha_s^2)/(besseli(0, alpha_s,1) + besseli(2, alpha_s,1)) + ((besselj(0, alpha_tilde_s,1) + besselj(2, alpha_tilde_s,1))*alpha_tilde_s^2)/(besselj(0, alpha_tilde_s,1) - besselj(2, alpha_tilde_s,1)))*W_TEprime_sf)/(Yf*Z0*alpha_0*alpha_tilde_sf*W_TE_sf) ;
            F4(ss,ss)          =     ((-1i/2)*(besseli(0, alpha_s,1) - besseli(2, alpha_s,1))*alpha_s*es*W_TM_sf)/((besseli(0, alpha_s,1) + besseli(2, alpha_s,1))*Z0*alpha_0*W_TE_sf) - ((1i/2)*(besseli(0, alpha_s,1) - besseli(2, alpha_s,1))*alpha_s*delta_s0*es*W_TM_sf)/((besseli(0, alpha_s,1) + besseli(2, alpha_s,1))*Z0*alpha_0*W_TE_sf) + (1i*sqrt(2)*((besseli(0, alpha_s,1) - besseli(2, alpha_s,1))/(2*(besseli(0, alpha_s,1) + besseli(2, alpha_s,1))) - (besselj(0, alpha_tilde_s,1) + besselj(2, alpha_tilde_s,1))/(2*(besselj(0, alpha_tilde_s,1) - besselj(2, alpha_tilde_s,1))))*alpha_f^2*alpha_s*sqrt(es)*W_TM_sf)/(Z0*alpha_0*alpha_tilde_sf^2*W_TE_sf) + (1i*(besselj(0, alpha_tilde_s,1) + besselj(2, alpha_tilde_s,1))*alpha_0*alpha_s*sqrt(es)*W_TM_sf)/(sqrt(2)*(besselj(0, alpha_tilde_s,1) - besselj(2, alpha_tilde_s,1))*Z0*alpha_tilde_sf^2*W_TE_sf) - ((1i/2)*(besseli(0, alpha_s,1) - besseli(2, alpha_s,1))*alpha_s^3*es*W_TM_sf)/((besseli(0, alpha_s,1) + besseli(2, alpha_s,1))*Z0*alpha_0*alpha_tilde_sf^2*W_TE_sf);
            
            G1(ss,pp)          =     1/(sqrt(-1 + beta_p^2)*beta_ps^2);
            G2(ss,pp)          =     1/(beta_p^2*sqrt(-1 + beta_p^2)*beta_ps^2*(alpha_0^2 - beta_ps^2));
            G3(ss,pp)          =     alpha_tilde_p/(alpha_p^2*(alpha_0^2 - alpha_ps^2));
            G4(ss,pp)          =     1/(sqrt(-1 + beta_p^2)*beta_ps^2);
            G5(ss,pp)          =     1/(sqrt(-1 + beta_p^2)*beta_ps^2*(alpha_0^2 - beta_ps^2));
             
            M1(pp,ss)          =     (2*(1i - ((2*1i)*(-1)^s)/(exp((1i*L*beta_tilde_p)/b)*(1 - exp(((-2*1i)*L*beta_tilde_p)/b))) + cot((L*beta_tilde_p)/b))*alpha_f*alpha_s*beta_p^2*(-(beta_ps^2*es*W_TEprime_sf) + alpha_0^2*(2*(-1 + delta_s0)*W_TEprime_sf + es*W_TEprime_sf)))/(beta_ps^2*(alpha_0^2 - beta_ps^2)*alpha_tilde_sf);
            M2(pp,ss)          =     ((1 - (-1)^s/exp((1i*L*alpha_tilde_p)/b))*sqrt(L*es)*W_TM_sf)/(alpha_0^2 - alpha_ps^2);
            M3(pp,ss)          =     (1i*sqrt(2)*(1i - ((2*1i)*(-1)^s)/(exp((1i*L*beta_tilde_p)/b)*(1 - exp(((-2*1i)*L*beta_tilde_p)/b))) + cot((L*beta_tilde_p)/b))*Yf*sqrt(es)*(2*alpha_f^2*alpha_0^2*alpha_s^2*(-1 + delta_s0) + beta_p^2*beta_ps^2*es*(alpha_s^2 + (1 + delta_s0)*alpha_tilde_sf^2) - alpha_0^2*(alpha_s^2*(2*beta_ps^2*(-1 + delta_s0) + beta_p^2*es) + beta_p^2*(1 + delta_s0)*es*alpha_tilde_sf^2))*W_TM_sf)/(beta_ps^2*(alpha_0^2 - beta_ps^2)*alpha_tilde_sf^2);
            M4(pp,ss)          =     ((2*1i)*(1i*(-1)^s - (2*1i)/(exp((1i*L*beta_tilde_p)/b)*(1 - exp(((-2*1i)*L*beta_tilde_p)/b))) + (-1)^s*cot((L*beta_tilde_p)/b))*alpha_f*alpha_s*beta_p^2*(-(beta_ps^2*es*W_TEprime_sf) + alpha_0^2*(2*(-1 + delta_s0)*W_TEprime_sf + es*W_TEprime_sf)))/(beta_ps^2*(alpha_0^2 - beta_ps^2)*alpha_tilde_sf);
            M5(pp,ss)          =     (sqrt(2)*(1i*(-1)^s - (2*1i)/(exp((1i*L*beta_tilde_p)/b)*(1 - exp(((-2*1i)*L*beta_tilde_p)/b))) + (-1)^s*cot((L*beta_tilde_p)/b))*Yf*sqrt(es)*(-2*alpha_f^2*alpha_0^2*alpha_s^2*(-1 + delta_s0) - beta_p^2*beta_ps^2*es*(alpha_s^2 + (1 + delta_s0)*alpha_tilde_sf^2) + alpha_0^2*(alpha_s^2*(2*beta_ps^2*(-1 + delta_s0) + beta_p^2*es) + beta_p^2*(1 + delta_s0)*es*alpha_tilde_sf^2))*W_TM_sf)/(beta_ps^2*(alpha_0^2 - beta_ps^2)*alpha_tilde_sf^2);
            
              
            N1(pp,pp)          =     ((-1i/4)*sqrt(pi)*beta_p^2)/(sqrt(L)*(1i + cot((L*beta_tilde_p)/b))*Yf*Z0*alpha_0*sqrt(-1 + beta_p^2)*beta_tilde_p);
            N2(pp,pp)          =     (sqrt(pi/2)*alpha_p^2)/L;
            N3(pp,pp)          =     (sqrt(pi)*beta_p^2)/(4*sqrt(L)*(1i + cot((L*beta_tilde_p)/b))*Yf*Z0*alpha_0*sqrt(-1 + beta_p^2)*beta_tilde_p);
            
            T1(ss,ss)          =     b*sqrt(2/pi);
            T2(ss,ss)          =     b*sqrt(2/pi)*alpha_0^2*alpha_s^2;
            T3(ss,ss)          =     (b*sqrt(2/pi)*alpha_0)/Z0;  
            T4(ss,ss)          =     (b*sqrt(2/pi)*alpha_s)/(sqrt(L/es)*W_TE_sf);
            T5(ss,ss)          =     (2*b*alpha_0^2*alpha_s)/(sqrt(L)*sqrt(pi)*W_TE_sf);
            
            % Matrices for impedance
            
            Z_integral_TM_sx(1,pp)        =    (1i*b*(beta*alpha_0 + alpha_tilde_p))/(sqrt(2*pi)*besselj(0, alpha_p)*alpha_p^2*(alpha_b + alpha_tilde_p));
            Z_integral_TE_sx(1,pp)        =    ((-1i)*b*(alpha_0 + beta*beta_tilde_p))/(sqrt(2*pi)*besselj(0, beta_p)*Y0*beta_p^2*sqrt(-1 + beta_p^2)*(alpha_b + beta_tilde_p));
            Z_integral_TM_dx(1,pp)        =    (1i*b*exp((1i*L*alpha_b)/b)*(-(beta*alpha_0) + alpha_tilde_p))/(sqrt(2*pi)*besselj(0, alpha_p)*alpha_p^2*(alpha_b - alpha_tilde_p));
            Z_integral_TE_dx(1,pp)        =    (1i*b*exp((1i*L*alpha_b)/b)*(alpha_0 - beta*beta_tilde_p))/(sqrt(2*pi)*besselj(0, beta_p)*Y0*beta_p^2*sqrt(-1 + beta_p^2)*(alpha_b - beta_tilde_p));
            
            Z_CTMsol(1,pp)                =    (b*((((2*1i)*exp((1i*L*alpha_b)/b - (1i*L*alpha_tilde_p)/b))/(1 - exp(((-2*1i)*L*alpha_tilde_p)/b)) - cot((L*alpha_tilde_p)/b))*alpha_b*alpha_p^2 + 1i*alpha_b^2*alpha_tilde_p - 1i*alpha_0^2*alpha_tilde_p))/(sqrt(2*pi)*besselj(2, alpha_p)*alpha_b*alpha_p^2*(alpha_b^2 - alpha_0^2 + alpha_p^2));
            Z_DTMsol(1,pp)                =    (b*(((-2*1i*exp((-1i*L*alpha_tilde_p)/b))/((1 - exp(((-2*1i)*L*alpha_tilde_p)/b))) + exp((1i*L*alpha_b)/b)*cot((L*alpha_tilde_p)/b))*alpha_b*alpha_p^2 + 1i*exp((1i*L*alpha_b)/b)*alpha_b^2*alpha_tilde_p - 1i*exp((1i*L*alpha_b)/b)*alpha_0^2*alpha_tilde_p))/(sqrt(2*pi)*besselj(2, alpha_p)*alpha_b*alpha_p^2*(alpha_b^2 - alpha_0^2 + alpha_p^2));
            Z_ATM1sol(1,ss)               =    (1i*b*(-1 + (-1)^s*exp((1i*L*alpha_b)/b))*(-alpha_0^2 + alpha_s^2)*(alpha_tilde_s^(-2) - 1/(2*besselj(1, alpha_tilde_s)*alpha_tilde_s))*W_TM_sf)/((alpha_b^2 - alpha_s^2)*sqrt(L/es));
            
            Z_CTEsol(1,pp)                =    (b*sqrt(L)*alpha_0*(1i*alpha_b + (-exp((1i*L*alpha_b)/b)*csch((L*beta_p)/b) + coth((L*beta_p)/b))*beta_p))/(sqrt(2*pi)*besselj(0, beta_p)*Y0*beta_p^2*sqrt(L*(-1 + beta_p^2))*(alpha_b^2 + beta_p^2));
            Z_DTEsol(1,pp)                =    -((b*alpha_0*(L*csch((L*beta_p)/b)*beta_p + 1i*exp((1i*L*alpha_b)/b)*L*(alpha_b + 1i*coth((L*beta_p)/b)*beta_p)))/(sqrt(L)*sqrt(2*pi)*besselj(0, beta_p)*Y0*beta_p^2*sqrt(L*(-1 + beta_p^2))*(alpha_b^2 + beta_p^2)));
            Z_ATEsol(1,ss)                =    exp(-abs(real(alpha_s)))*(sqrt(2)*b*(-1 + (-1)^s*exp((1i*L*alpha_b)/b))*alpha_f*alpha_s*(-1 + delta_s0)*W_TEprime_sf)/(sqrt(L)*(besseli(0, alpha_s,1) + besseli(2, alpha_s,1))*Yf*(-alpha_b^2 + alpha_s^2)*alpha_tilde_sf);
            Z_ATM2sol(1,ss)               =    ((-1i)*b*(-1 + (-1)^s*exp((1i*L*alpha_b)/b))*((-1 + exp(-abs(real(alpha_s)))*(besseli(0, alpha_s,1) + besseli(2, alpha_s,1))^(-1))*alpha_f^2 + alpha_s^2)*sqrt(es)*W_TM_sf)/(sqrt(L)*(alpha_b^2 - alpha_s^2)*alpha_tilde_sf^2);
            
            Z_CTEirr(1,pp)                =    (b*Z0*alpha_0*((-coth((L*beta_p)/b) + exp((1i*L*alpha_b)/b)*csch((L*beta_p)/b))*alpha_b + 1i*beta_p))/(sqrt(2*pi)*besselj(0, beta_p)*alpha_b*beta_p*sqrt(-1 + beta_p^2)*(alpha_b^2 + beta_p^2));
            Z_DTEirr(1,pp)                =    (b*beta*Z0*(csch((L*beta_p)/b)*alpha_b - exp((1i*L*alpha_b)/b)*(coth((L*beta_p)/b)*alpha_b + 1i*beta_p)))/(sqrt(2*pi)*besselj(0, beta_p)*beta_p*sqrt(-1 + beta_p^2)*(alpha_b^2 + beta_p^2));
            Z_ATMirr(1,ss)                =    exp(-abs(real(alpha_s)))*((-1i/2)*b*(-1 + (-1)^s*exp((1i*L*alpha_b)/b))*beta*alpha_b*es^(3/2)*(alpha_s^2 + (1 + delta_s0)*alpha_tilde_sf^2)*W_TM_sf)/(sqrt(L)*(-2*besseli(0, alpha_s,1) + (besseli(0, alpha_s,1) - besseli(2, alpha_s,1)))*alpha_0*(alpha_b^2 - alpha_s^2)*alpha_tilde_sf^2);
            Z_ATEirr(1,ss)                =    -exp(-abs(real(alpha_s)))*((b*(-1 + (-1)^s*exp((1i*L*alpha_0)/(b*beta)))*beta^2*alpha_f*alpha_s*sqrt(es)*W_TEprime_sf)/(sqrt(2)*(-besseli(0, alpha_s,1) - besseli(2, alpha_s,1))*Yf*(-alpha_0^2 + beta^2*alpha_s^2)*sqrt(L/es)*alpha_tilde_sf));
            
        end

    end
    
    % Accuracy check
    
    if ~isfinite(exp((1i*L*beta_tilde_p)/b)); warning('Infinity encountered, results may not be accurate....'); end

    % Inversions

    FF1=F4 - T4*G4*N1*M3 + T5*G5*N1*M3 + IIss*T4*G4*N3*M5 - IIss*T5*G5*N3*M5;
    FF2=F3 - T4*G4*N1*M1 + T5*G5*N1*M1 + IIss*T4*G4*N3*M4 - IIss*T5*G5*N3*M4;
    UU1=FF2\FF1;
    UU2=F2 - F1*UU1 - T1*G1*N1*M3 - T2*G2*N1*M3 + T3*G3*N2*M2 + IIss*T1*G1*N3*M5 + IIss*T2*G2*N3*M5 + T1*G1*N1*M1*UU1 + T2*G2*N1*M1*UU1 - IIss*T1*G1*N3*M4*UU1 - IIss*T2*G2*N3*M4*UU1 + IIss*T3*G3*N2*M2*IIss;
        
    A_TM=-UU2\(B);
    A_TE=-UU1*(A_TM);
    C_TM=N2*M2*(A_TM);
    D_TM=N2*M2*IIss*(A_TM);
    C_TE=N1*M1*(A_TE) + N1*M3*(A_TM);
    D_TE=N3*M4*(A_TE) + N3*M5*(A_TM);

    
    Zsx=Z_integral_TM_sx*(C_TM)+Z_integral_TE_sx*(C_TE);
    Zdx=Z_integral_TM_dx*(D_TM)+Z_integral_TE_dx*(D_TE);
    Zcav_sol=Z_CTMsol*(C_TM)+Z_DTMsol*(D_TM)+Z_CTEsol*(C_TE)+Z_DTEsol*(D_TE)+(Z_ATM1sol+Z_ATM2sol)*(A_TM)+Z_ATEsol*(A_TE);
    Zcav_irr=Z_CTEirr*(C_TE)+Z_DTEirr*(D_TE)+Z_ATMirr*(A_TM)+Z_ATEirr*(A_TE);

    Z(ff)=Zsx+Zdx+Zcav_sol+Zcav_irr;
 

            
    delta_skin=sqrt(2/(2*pi*f*muf*sigma));
    Zrew(ff)=sl/(2*pi*f)*beta*(1+1i)/(delta_skin*sigma*b^3*pi);

    if ismember(f,f_post)
        A_TM_post=[A_TM_post,A_TM];
        C_TM_post=[C_TM_post,C_TM];
        D_TM_post=[D_TM_post,D_TM];
        A_TE_post=[A_TE_post,A_TE];
        C_TE_post=[C_TE_post,C_TE];
        D_TE_post=[D_TE_post,D_TE];       
    end
end


fco=Jzeroprime(1)*sl/2/pi/b/1e9;
TEmodes_WG=Jzeroprime*sl/2/pi/b/1e9;
TMmodes_WG=Jzero*sl/2/pi/b/1e9;
TEmodes_C=[];
for s=0:S-1
    if s~=0
        TEmodes_C=[TEmodes_C,sl/2/pi/c/1e9*sqrt(Jzeroprime.^2+(s*pi/L*c).^2)];
    end
end

TMmodes_C=[];
for s=0:S-1
    TMmodes_C=[TMmodes_C,sl/2/pi/c/1e9*sqrt(Jzero.^2+(s*pi/L*c).^2)];
end

% disp(TMmodes_C)
% disp(TEmodes_C)
% disp(TMmodes_WG);
% disp(TEmodes_WG);
%
% set(0,'defaultlinelinewidth',2)
% figure(2)
% name=Z;
% loglog(f_vec,real(name)/L,'-r'); hold on;
% loglog(f_vec,imag(name)/L,'--g'); 
% plot((f_vec),real(Zrew),'-k','markerfacecolor','k');
% plot((f_vec),imag(Zrew),'--r');
% hold off;
% legend('MMM Re(Z_{dip})','MMM Im(Z_{dip})','ReW Re(Z_{dip})','ReW Im(Z_{dip})');
% grid on;
% xlabel('f');
% ylabel('Ohm /m');
% title(['f_c_o =',num2str(fco),'GHz'])
% vline(fco*1e9,'--r',{'f_c_o'})
% xlim([1 1e10])
% ylim([1 1e10])
    
%%
Fields_vec=[];

for f=f_post
    if fields==1
    disp(['Postprocessing Fields  at ',num2str(f),'...']);
    
    % Material properties
    
    omega=2*pi*f;
    epsilonf=epsilon0*(eval(e_r)-1i*sigma/(omega*epsilon0)); % Relative dielectric 
    muf=mu0*eval(mu_r);
    Zf=sqrt(muf/epsilonf);
    Yf=1/Zf;
   

    % Kappas
    k0=omega/sl;
    kf=omega*sqrt(muf*epsilonf);
    kb=omega/(beta*sl);
  
    % Meshing
    
    Dzz=L/5;
    Drr=b/5;
   
    PEz=[];
    PEr=[];
    PEphi=[];
    PHz=[];
    PHr=[];
    PHphi=[];

    rr_mesh=[Drr/2:Drr:c];
    zz_mesh=[-L-Dzz/2:Dzz:2*L];
    
    % Coefficients
    
    [~,indf]=ismember(f,f_post);
    A_TM=A_TM_post(:,indf);
    C_TM=C_TM_post(:,indf);
    D_TM=D_TM_post(:,indf);
    A_TE=A_TE_post(:,indf);
    C_TE=C_TE_post(:,indf);
    D_TE=D_TE_post(:,indf);    
    
    for r_ind=1:length(rr_mesh)
        rr=rr_mesh(r_ind);
        disp([num2str(rr),'/',num2str(c)]);
        
       for z_ind=1:length(zz_mesh)
           zz=zz_mesh(z_ind);
                 
                Ez=0;
                Er=0;
                Ephi=0;
                Hz=0;
                Hr=0;
                Hphi=0;
                Ez_b=0;
                Er_b=0;
                Ephi_b=0;
                Hz_b=0;
                Hr_b=0;
                Hphi_b=0;
                
                % LOAD
                
                if rr>b && rr<=c
                if zz>0 && zz<=L
                        
                    for s=0:S-1;

                        ss=s+1; %indici per riempire le matrici
                        if s==0 es=1;  else es=2;    end
                        if s==0 delta_s0=1; else delta_s0=0; end

                        alpha_0           =       k0*b;
                        alpha_02          =       alpha_0^2;
                        alpha_f           =       kf*b;
                        alpha_f2          =       alpha_f^2;
                        alpha_b           =       kb*b;
                        alpha_b2          =       alpha_b^2;
                        alpha_s           =       s*pi/(L/b);
                        alpha_s2          =       alpha_s^2;
                        alpha_tilde_s2    =       alpha_0^2-alpha_s^2;
                        alpha_tilde_sf2   =       alpha_f^2-alpha_s^2;
                        alpha_tilde_sf    =       sqrt(alpha_tilde_sf2);
                        alpha_tilde_s     =       (sqrt(alpha_tilde_s2));

                        beta_s            =       alpha_s;
                        beta_s2           =       beta_s^2;
                        beta_tilde_s2     =       alpha_02-beta_s2;
                        beta_tilde_s      =       (sqrt(beta_tilde_s2));

                        if sigma==0
                            alpha_tilde_sf    =       conj(alpha_tilde_sf);
                        end

                        in                 =        alpha_tilde_sf;
                        W_TE_sf            =        (exp((1i*(b - rr)*alpha_tilde_sf)/b)*(besselh(0,1,(c*alpha_tilde_sf)/b,1) - besselh(2,1,(c*alpha_tilde_sf)/b,1))*besselh(1,2,(rr*alpha_tilde_sf)/b,1) - exp((1i*(b - 2*c + rr)*alpha_tilde_sf)/b)*besselh(1,1,(rr*alpha_tilde_sf)/b,1)*(besselh(0,2,(c*alpha_tilde_sf)/b,1) - besselh(2,2,(c*alpha_tilde_sf)/b,1)))/((besselh(0,1,(c*alpha_tilde_sf)/b,1) - besselh(2,1,(c*alpha_tilde_sf)/b,1))*besselh(1,2,alpha_tilde_sf,1));
                        W_TEprime_sf       =        (-(exp((1i*(b - 2*c + rr)*alpha_tilde_sf)/b)*besselh(0,1,(rr*alpha_tilde_sf)/b,1)*(besselh(0,2,(c*alpha_tilde_sf)/b,1) - besselh(2,2,(c*alpha_tilde_sf)/b,1))) + exp((1i*(b - 2*c + rr)*alpha_tilde_sf)/b)*besselh(2,1,(rr*alpha_tilde_sf)/b,1)*(besselh(0,2,(c*alpha_tilde_sf)/b,1) - besselh(2,2,(c*alpha_tilde_sf)/b,1)) + exp((1i*(b - rr)*alpha_tilde_sf)/b)*(besselh(0,1,(c*alpha_tilde_sf)/b,1) - besselh(2,1,(c*alpha_tilde_sf)/b,1))*(besselh(0,2,(rr*alpha_tilde_sf)/b,1) - besselh(2,2,(rr*alpha_tilde_sf)/b,1)))/(2*(besselh(0,1,(c*alpha_tilde_sf)/b,1) - besselh(2,1,(c*alpha_tilde_sf)/b,1))*besselh(1,2,alpha_tilde_sf,1));
                        W_TM_sf            =        (-(exp((1i*(b - 2*c + rr)*alpha_tilde_sf)/b)*besselh(1,1,(rr*alpha_tilde_sf)/b,1)*besselh(1,2,(c*alpha_tilde_sf)/b,1)) + exp((1i*(b - rr)*alpha_tilde_sf)/b)*besselh(1,1,(c*alpha_tilde_sf)/b,1)*besselh(1,2,(rr*alpha_tilde_sf)/b,1))/(besselh(1,1,(c*alpha_tilde_sf)/b,1)*besselh(1,2,alpha_tilde_sf,1));
                        W_TMprime_sf       =        (-(exp((1i*(b - 2*c + rr)*alpha_tilde_sf)/b)*(besselh(0,1,(rr*alpha_tilde_sf)/b,1) - besselh(2,1,(rr*alpha_tilde_sf)/b,1))*besselh(1,2,(c*alpha_tilde_sf)/b,1)) + exp((1i*(b - rr)*alpha_tilde_sf)/b)*besselh(1,1,(c*alpha_tilde_sf)/b,1)*(besselh(0,2,(rr*alpha_tilde_sf)/b,1) - besselh(2,2,(rr*alpha_tilde_sf)/b,1)))/(2*besselh(1,1,(c*alpha_tilde_sf)/b,1)*besselh(1,2,alpha_tilde_sf,1));
            
                        Er=Er +         -(cos(ph)*sin((zz*alpha_s)/b)*A_TM(ss)*alpha_s*W_TMprime_sf/(sqrt(L/es)*alpha_tilde_sf))+(-1i)*sqrt(2)*b*cos(ph)*sin((zz*alpha_s)/b)*A_TE(ss)*alpha_f*W_TE_sf/(sqrt(L)*rr*Yf*alpha_tilde_sf^2);
                        Ephi=Ephi +     b*sin(ph)*sin((zz*alpha_s)/b)*A_TM(ss)*alpha_s*W_TM_sf/(rr*sqrt(L/es)*alpha_tilde_sf^2)+1i*sqrt(2)*sin(ph)*sin((zz*alpha_s)/b)*A_TE(ss)*alpha_f*W_TEprime_sf/(sqrt(L)*Yf*alpha_tilde_sf);
                        Ez=Ez +         cos(ph)*cos((zz*alpha_s)/b)*A_TM(ss)*W_TM_sf/sqrt(L/es)+0;
                        Hr=Hr +         (-1i)*b*cos((zz*alpha_s)/b)*sin(ph)*A_TM(ss)*alpha_f*W_TM_sf/(rr*Zf*sqrt(L/es)*alpha_tilde_sf^2)+sqrt(2)*cos((zz*alpha_s)/b)*sin(ph)*A_TE(ss)*alpha_s*W_TEprime_sf/(sqrt(L)*alpha_tilde_sf);  
                        Hphi=Hphi +     (-1i)*cos(ph)*cos((zz*alpha_s)/b)*A_TM(ss)*alpha_f*W_TMprime_sf/(Zf*sqrt(L/es)*alpha_tilde_sf)+sqrt(2)*b*cos(ph)*cos((zz*alpha_s)/b)*A_TE(ss)*alpha_s*W_TE_sf/(sqrt(L)*rr*alpha_tilde_sf^2); 
                        Hz=Hz +         0+sqrt(2)*sin(ph)*sin((zz*alpha_s)/b)*A_TE(ss)*W_TE_sf/sqrt(L); 

                        Ez_b=0;
                        Er_b=0;
                        Ephi_b=0;
                        Hz_b=0;
                        Hr_b=0;
                        Hphi_b=0;
          
                 
                    end
                                        
                end
                end
                    
                % PIPE SX
                    
                if rr<b && zz<0
                
                    
                    for p=1:P
                        pp=p;
                        alpha_0           =       k0*b;
                        alpha_02          =       alpha_0^2;
                        alpha_b           =       kb*b;
                        alpha_b2          =       alpha_b^2;
                        alpha_p           =       Jzero(p);
                        alpha_p2          =       alpha_p^2;
                        alpha_tilde_p2    =       alpha_0^2-alpha_p^2;
                        alpha_tilde_p     =       conj(sqrt(alpha_tilde_p2));

                        beta_p            =       Jzeroprime(p);
                        beta_p2           =       beta_p^2;
                        beta_tilde_p2     =       alpha_02-beta_p2;
                        beta_tilde_p      =       conj(sqrt(beta_tilde_p2));
                        
                        Er=Er      +(1i*exp((1i*zz*alpha_tilde_p)/b)*(besselj(0, (rr*alpha_p)/b) - besselj(2, (rr*alpha_p)/b))*cos(ph)*C_TM(pp)*alpha_tilde_p)/(sqrt(2*pi)*besselj(0, alpha_p)*alpha_p^2)  +...
                                    ((-1i)*b*exp((1i*zz*beta_tilde_p)/b)*sqrt(2/pi)*besselj(1, (rr*beta_p)/b)*cos(ph)*C_TE(pp)*alpha_0)/(rr*besselj(0, beta_p)*Y0*beta_p^3*sqrt(-1 + beta_p^2));       
                        Ephi=Ephi  +((-1i)*b*exp((1i*zz*alpha_tilde_p)/b)*sqrt(2/pi)*besselj(1, (rr*alpha_p)/b)*sin(ph)*C_TM(pp)*alpha_tilde_p)/(rr*besselj(0, alpha_p)*alpha_p^3) +    ...
                                    ((-1i)*b*exp((1i*zz*beta_tilde_p)/b)*sqrt(2/pi)*C_TE(pp)*alpha_0*(-(besselj(0, (rr*beta_p)/b)*sin(ph)*beta_p)/(2*b) + (besselj(2, (rr*beta_p)/b)*sin(ph)*beta_p)/(2*b)))/(besselj(0, beta_p)*Y0*beta_p^3*sqrt(-1 + beta_p^2));
                                
                        Ez=Ez      +(exp((1i*zz*alpha_tilde_p)/b)*sqrt(2/pi)*besselj(1, (rr*alpha_p)/b)*cos(ph)*C_TM(pp))/(besselj(0, alpha_p)*alpha_p)   +...
                                    0;
                        Hr=Hr      +((-1i)*b*exp((1i*zz*alpha_tilde_p)/b)*sqrt(2/pi)*besselj(1, (rr*alpha_p)/b)*sin(ph)*C_TM(pp)*alpha_0)/(rr*besselj(0, alpha_p)*Z0*alpha_p^3)  +...
                                    (1i*exp((1i*zz*beta_tilde_p)/b)*(besselj(0, (rr*beta_p)/b) - besselj(2, (rr*beta_p)/b))*sin(ph)*C_TE(pp)*beta_tilde_p)/(sqrt(2*pi)*besselj(0, beta_p)*beta_p^2*sqrt(-1 + beta_p^2));
                        Hphi=Hphi  +(1i*b*exp((1i*zz*alpha_tilde_p)/b)*sqrt(2/pi)*C_TM(pp)*alpha_0*(-(besselj(0, (rr*alpha_p)/b)*cos(ph)*alpha_p)/(2*b) + (besselj(2, (rr*alpha_p)/b)*cos(ph)*alpha_p)/(2*b)))/(besselj(0, alpha_p)*Z0*alpha_p^3) +...
                                    (1i*b*exp((1i*zz*beta_tilde_p)/b)*sqrt(2/pi)*besselj(1, (rr*beta_p)/b)*cos(ph)*C_TE(pp)*beta_tilde_p)/(rr*besselj(0, beta_p)*beta_p^3*sqrt(-1 + beta_p^2));
                        Hz=Hz      +(exp((1i*zz*beta_tilde_p)/b)*sqrt(2/pi)*besselj(1, (rr*beta_p)/b)*sin(ph)*C_TE(pp))/(besselj(0, beta_p)*beta_p*sqrt(-1 + beta_p^2))+...
                                    0;
                                
                        Ez_b=0;
                        Er_b=0;
                        Ephi_b=0;
                        Hz_b=0;
                        Hr_b=0;
                        Hphi_b=0;
                    end
               
                end
                
                % PIPE DX
                
                if rr<b && zz>L
               
                    for p=1:P
                        pp=p;
                        alpha_0           =       k0*b;
                        alpha_02          =       alpha_0^2;
                        alpha_b           =       kb*b;
                        alpha_b2          =       alpha_b^2;
                        alpha_p           =       Jzero(p);
                        alpha_p2          =       alpha_p^2;
                        alpha_tilde_p2    =       alpha_0^2-alpha_p^2;
                        alpha_tilde_p     =       conj(sqrt(alpha_tilde_p2));

                        beta_p            =       Jzeroprime(p);
                        beta_p2           =       beta_p^2;
                        beta_tilde_p2     =       alpha_02-beta_p2;
                        beta_tilde_p      =       conj(sqrt(beta_tilde_p2));

                        Er=Er     +((-1i)*(besselj(0, (rr*alpha_p)/b) - besselj(2, (rr*alpha_p)/b))*cos(ph)*D_TM(pp)*alpha_tilde_p)/(exp((1i*(-L + zz)*alpha_tilde_p)/b)*sqrt(2*pi)*besselj(0, alpha_p)*alpha_p^2)+...
                                   ((-1i)*b*sqrt(2/pi)*besselj(1, (rr*beta_p)/b)*cos(ph)*D_TE(pp)*alpha_0)/(exp((1i*(-L + zz)*beta_tilde_p)/b)*rr*besselj(0, beta_p)*Y0*beta_p^3*sqrt(-1 + beta_p^2));
                        Ephi=Ephi +(1i*b*sqrt(2/pi)*besselj(1, (rr*alpha_p)/b)*sin(ph)*D_TM(pp)*alpha_tilde_p)/(exp((1i*(-L + zz)*alpha_tilde_p)/b)*rr*besselj(0, alpha_p)*alpha_p^3)+...
                                   ((-1i)*b*sqrt(2/pi)*D_TE(pp)*alpha_0*(-(besselj(0, (rr*beta_p)/b)*sin(ph)*beta_p)/(2*b) + (besselj(2, (rr*beta_p)/b)*sin(ph)*beta_p)/(2*b)))/(exp((1i*(-L + zz)*beta_tilde_p)/b)*besselj(0, beta_p)*Y0*beta_p^3*sqrt(-1 + beta_p^2));
                        Ez=Ez     +(sqrt(2/pi)*besselj(1, (rr*alpha_p)/b)*cos(ph)*D_TM(pp))/(exp((1i*(-L + zz)*alpha_tilde_p)/b)*besselj(0, alpha_p)*alpha_p)+...
                                   0*D_TE(pp);
                        Hr=Hr     +((-1i)*b*sqrt(2/pi)*besselj(1, (rr*alpha_p)/b)*sin(ph)*D_TM(pp)*alpha_0)/(exp((1i*(-L + zz)*alpha_tilde_p)/b)*rr*besselj(0, alpha_p)*Z0*alpha_p^3)+...
                                   ((-1i)*(besselj(0, (rr*beta_p)/b) - besselj(2, (rr*beta_p)/b))*sin(ph)*D_TE(pp)*beta_tilde_p)/(exp((1i*(-L + zz)*beta_tilde_p)/b)*sqrt(2*pi)*besselj(0, beta_p)*beta_p^2*sqrt(-1 + beta_p^2));
                        Hphi=Hphi +(1i*b*sqrt(2/pi)*D_TM(pp)*alpha_0*(-(besselj(0, (rr*alpha_p)/b)*cos(ph)*alpha_p)/(2*b) + (besselj(2, (rr*alpha_p)/b)*cos(ph)*alpha_p)/(2*b)))/(exp((1i*(-L + zz)*alpha_tilde_p)/b)*besselj(0, alpha_p)*Z0*alpha_p^3)+...
                                   ((-1i)*b*sqrt(2/pi)*besselj(1, (rr*beta_p)/b)*cos(ph)*D_TE(pp)*beta_tilde_p)/(exp((1i*(-L + zz)*beta_tilde_p)/b)*rr*besselj(0, beta_p)*beta_p^3*sqrt(-1 + beta_p^2)); %ok%
                        Hz=Hz     +(sqrt(2/pi)*besselj(1, (rr*beta_p)/b)*sin(ph)*D_TE(pp))/(exp((1i*(-L + zz)*beta_tilde_p)/b)*besselj(0, beta_p)*beta_p*sqrt(-1 + beta_p^2))+...
                                   0*D_TM(pp); %ok%
                               
                        Ez_b=0;
                        Er_b=0;
                        Ephi_b=0;
                        Hz_b=0;
                        Hr_b=0;
                        Hphi_b=0;
                    end
               
                end
                 
                % CAVITY
                
                if rr<b && zz>0 && zz<L
                       
                    for s=0:S-1
                        ss=s+1;
                        if s==0 es=1;  else es=2;    end
                        if s==0 delta_s0=1; else delta_s0=0; end
                        for p=1:P
                            pp=p;
                        alpha_0           =       k0*b;
                        alpha_02          =       alpha_0^2;
                        alpha_f           =       kf*b;
                        alpha_f2          =       alpha_f^2;
                        alpha_b           =       kb*b;
                        alpha_b2          =       alpha_b^2;
                        alpha_p           =       Jzero(p);
                        alpha_s           =       s*pi/(L/b);
                        alpha_s2          =       alpha_s^2;
                        alpha_p2          =       alpha_p^2;
                        alpha_ps          =       sqrt(alpha_p^2+alpha_s^2);
                        alpha_ps2         =       alpha_ps^2;
                        alpha_tilde_s2    =       alpha_0^2-alpha_s^2;
                        alpha_tilde_p2    =       alpha_0^2-alpha_p^2;
                        alpha_tilde_sf2   =       alpha_f^2-alpha_s^2;
                        alpha_tilde_sf    =       sqrt(alpha_tilde_sf2);
                        alpha_tilde_p     =       conj(sqrt(alpha_tilde_p2));
                        alpha_tilde_s     =       (sqrt(alpha_tilde_s2));

                        beta_s            =       alpha_s;
                        beta_s2           =       beta_s^2;
                        beta_p            =       Jzeroprime(p);
                        beta_p2           =       beta_p^2;
                        beta_ps2          =       beta_p2+beta_s2;
                        beta_tilde_p2     =       alpha_02-beta_p2;
                        beta_tilde_s2     =       alpha_02-beta_s2;
                        beta_tilde_p      =       conj(sqrt(beta_tilde_p2));
                        beta_tilde_s      =       (sqrt(beta_tilde_s2));
                        beta_ps           =       sqrt(beta_ps2);

                        if sigma==0
                            alpha_tilde_sf    =       conj(alpha_tilde_sf);
                        end

                        in                 =        alpha_tilde_sf;
 
                        W_TE_sf            =        ((besselh(0,1,(c*alpha_tilde_sf)/b,1) - besselh(2,1,(c*alpha_tilde_sf)/b,1))*besselh(1,2,alpha_tilde_sf,1) - exp(((2*1i)*(b - c)*alpha_tilde_sf)/b)*besselh(1,1,alpha_tilde_sf,1)*(besselh(0,2,(c*alpha_tilde_sf)/b,1) - besselh(2,2,(c*alpha_tilde_sf)/b,1)))/((besselh(0,1,(c*alpha_tilde_sf)/b,1) - besselh(2,1,(c*alpha_tilde_sf)/b,1))*besselh(1,2,alpha_tilde_sf,1));
                        W_TEprime_sf       =        (besselh(0,1,(c*alpha_tilde_sf)/b,1)*(besselh(0,2,alpha_tilde_sf,1) - besselh(2,2,alpha_tilde_sf,1)) + besselh(2,1,(c*alpha_tilde_sf)/b,1)*(-besselh(0,2,alpha_tilde_sf,1) + besselh(2,2,alpha_tilde_sf,1)) - exp(((2*1i)*(b - c)*alpha_tilde_sf)/b)*(besselh(0,1,alpha_tilde_sf,1) - besselh(2,1,alpha_tilde_sf,1))*(besselh(0,2,(c*alpha_tilde_sf)/b,1) - besselh(2,2,(c*alpha_tilde_sf)/b,1)))/(2*(besselh(0,1,(c*alpha_tilde_sf)/b,1) - besselh(2,1,(c*alpha_tilde_sf)/b,1))*besselh(1,2,alpha_tilde_sf,1));
                        W_TM_sf            =        1 - (exp(((2*1i)*(b - c)*alpha_tilde_sf)/b)*besselh(1,1,alpha_tilde_sf,1)*besselh(1,2,(c*alpha_tilde_sf)/b,1))/(besselh(1,1,(c*alpha_tilde_sf)/b,1)*besselh(1,2,alpha_tilde_sf,1));
                        W_TMprime_sf       =        (-(((besselh(0,1,alpha_tilde_sf,1) - besselh(2,1,alpha_tilde_sf,1))*besselh(1,2,(c*alpha_tilde_sf)/b,1))/exp(((2*1i)*(-b + c)*alpha_tilde_sf)/b)) + besselh(1,1,(c*alpha_tilde_sf)/b,1)*(besselh(0,2,alpha_tilde_sf,1) - besselh(2,2,alpha_tilde_sf,1)))/(2*besselh(1,1,(c*alpha_tilde_sf)/b,1)*besselh(1,2,alpha_tilde_sf,1));
            
                        I_TM=(1i*b*alpha_0*((1i*b*C_TM(pp)*alpha_tilde_p)/(alpha_p^2*sqrt(L/es)) + (1i*(-1)^s*b*D_TM(pp)*alpha_tilde_p)/(alpha_p^2*sqrt(L/es)) + sqrt(2*pi)*A_TM(ss)*W_TM_sf))/(Z0*(alpha_0^2 - alpha_ps^2));
                        I_TE=(1i*b*alpha_0*((1i*sqrt(2)*b*C_TE(pp)*alpha_0*alpha_s)/(sqrt(L)*Y0*beta_p^2*beta_ps) - (1i*(-1)^s*sqrt(2)*b*D_TE(pp)*alpha_0*alpha_s)/(sqrt(L)*Y0*beta_p^2*beta_ps) - (sqrt(pi)*A_TM(ss)*alpha_s*(alpha_f^2 - beta_ps^2)*sqrt(es)*W_TM_sf)/(sqrt(-1 + beta_p^2)*beta_ps*alpha_tilde_sf^2) - (1i*sqrt(2*pi)*A_TE(ss)*alpha_f*beta_p^2*(-1 + delta_s0)*W_TEprime_sf)/(Yf*sqrt(-1 + beta_p^2)*beta_ps*alpha_tilde_sf)))/(Z0*(alpha_0^2 - beta_ps^2));
                        G_ps=(1i*b*Y0*((1i*b*C_TE(pp)*alpha_0)/(Y0*beta_p*beta_ps*sqrt(L/es)) - (1i*(-1)^s*b*D_TE(pp)*alpha_0)/(Y0*beta_p*beta_ps*sqrt(L/es)) - (sqrt(pi/2)*A_TM(ss)*beta_p*es*(alpha_s^2 + (1 + delta_s0)*alpha_tilde_sf^2)*W_TM_sf)/(sqrt(-1 + beta_p^2)*beta_ps*alpha_tilde_sf^2) - (1i*sqrt(pi)*A_TE(ss)*alpha_f*alpha_s*beta_p*sqrt(es)*W_TEprime_sf)/(Yf*sqrt(-1 + beta_p^2)*beta_ps*alpha_tilde_sf)))/alpha_0;
                        
                        V_TM=(((-1i)*Z0*alpha_ps)/alpha_0)*I_TM;
                        V_TE=(((-1i)*Z0*beta_ps)/alpha_0)*I_TE;
                        Er=Er +      V_TM*     -(((besselj(0, (rr*alpha_p)/b) - besselj(2, (rr*alpha_p)/b))*cos(ph)*sin((zz*alpha_s)/b)*alpha_s)/(b*sqrt(2*pi)*besselj(0, alpha_p)*alpha_ps*sqrt(L/es))) + V_TE*(2*besselj(1, (rr*beta_p)/b)*cos(ph)*sin((zz*alpha_s)/b))/(sqrt(pi)*rr*besselj(0, beta_p)*beta_p*sqrt(L*(-1 + beta_p^2)));
                        Ephi=Ephi +  V_TM*(sqrt(2/pi)*besselj(1, (rr*alpha_p)/b)*sin(ph)*sin((zz*alpha_s)/b)*alpha_s)/(rr*besselj(0, alpha_p)*alpha_p*alpha_ps*sqrt(L/es)) + V_TE*(2*sin((zz*alpha_s)/b)*(-(besselj(0, (rr*beta_p)/b)*sin(ph)*beta_p)/(2*b) + (besselj(2, (rr*beta_p)/b)*sin(ph)*beta_p)/(2*b)))/(sqrt(pi)*besselj(0, beta_p)*beta_p*sqrt(L*(-1 + beta_p^2)));
                        Ez=Ez +      V_TM* (sqrt(2/pi)*besselj(1, (rr*alpha_p)/b)*cos(ph)*cos((zz*alpha_s)/b)*alpha_p)/(b*besselj(0, alpha_p)*alpha_ps*sqrt(L/es)) + V_TE*0;
                        Hr=Hr +      I_TM* -((sqrt(2/pi)*besselj(1, (rr*alpha_p)/b)*cos((zz*alpha_s)/b)*sin(ph))/(rr*besselj(0, alpha_p)*alpha_p*sqrt(L/es)))+I_TE*  ((besselj(0, (rr*beta_p)/b) - besselj(2, (rr*beta_p)/b))*cos((zz*alpha_s)/b)*sin(ph)*alpha_s)/(b*sqrt(pi)*besselj(0, beta_p)*sqrt(L*(-1 + beta_p^2))*beta_ps)+G_ps*((besselj(0, (rr*beta_p)/b) - besselj(2, (rr*beta_p)/b))*cos((zz*alpha_s)/b)*sin(ph)*beta_p)/(b*sqrt(2*pi)*besselj(0, beta_p)*beta_ps*sqrt((L*(-1 + beta_p^2))/es));
                        Hphi=Hphi +  I_TM*(sqrt(2/pi)*cos((zz*alpha_s)/b)*(-(besselj(0, (rr*alpha_p)/b)*cos(ph)*alpha_p)/(2*b) + (besselj(2, (rr*alpha_p)/b)*cos(ph)*alpha_p)/(2*b)))/(besselj(0, alpha_p)*alpha_p*sqrt(L/es))+I_TE*(2*besselj(1, (rr*beta_p)/b)*cos(ph)*cos((zz*alpha_s)/b)*alpha_s)/(sqrt(pi)*rr*besselj(0, beta_p)*beta_p*sqrt(L*(-1 + beta_p^2))*beta_ps)+G_ps*(sqrt(2/pi)*besselj(1, (rr*beta_p)/b)*cos(ph)*cos((zz*alpha_s)/b))/(rr*besselj(0, beta_p)*beta_ps*sqrt((L*(-1 + beta_p^2))/es));
                        Hz=Hz +      I_TM*0+I_TE*(2*besselj(1, (rr*beta_p)/b)*sin(ph)*sin((zz*alpha_s)/b)*beta_p)/(b*sqrt(pi)*besselj(0, beta_p)*sqrt(L*(-1 + beta_p^2))*beta_ps)+G_ps*-((sqrt(2/pi)*besselj(1, (rr*beta_p)/b)*sin(ph)*sin((zz*alpha_s)/b)*alpha_s)/(b*besselj(0, beta_p)*beta_ps*sqrt((L*(-1 + beta_p^2))/es)));
                    
                        Ez_b=0;
                        Er_b=0;
                        Ephi_b=0;
                        Hz_b=0;
                        Hr_b=0;
                        Hphi_b=0;
                        end
                    end
                    
                    
                    
                
                end
                
                if rr>b
                    if zz>L || zz<0
                           Ez=nan;
                           Er=nan;
                           Ephi=nan;
                           Hz=nan;
                           Hr=nan;
                           Hphi=nan;

                    end
                end
                PEz(r_ind,z_ind)=Ez+Ez_b; % E Field vectors
                PEr(r_ind,z_ind)=Er+Er_b;
                PEphi(r_ind,z_ind)=Ephi+Ephi_b; 
                PHz(r_ind,z_ind)=Hz+Hz_b; % H Field vectors
                PHr(r_ind,z_ind)=Hr+Hr_b;
                PHphi(r_ind,z_ind)=Hphi+Hphi_b; 
       
       end % end z loop
    end % end r loop
      
 
    %% Plot Ez,r Hphi
   
    figure(); pcolor(zz_mesh,rr_mesh,(abs(PEr))); shading interp; hold on; contour(zz_mesh,rr_mesh,(abs(PEr))); hold off;
    xlabel('Length [m]'); ylabel('Radius [m]'); title(['Field E_r [V/m] @ f=',num2str(f)]); caxis([-max(max(abs(PEr))) max(max(abs(PEr)))]); colorbar;
    figure(); pcolor(zz_mesh,rr_mesh,(abs(PEz))); shading interp; hold on; contour(zz_mesh,rr_mesh,(abs(PEz))); hold off;
    xlabel('Length [m]'); ylabel('Radius [m]'); title(['Field E_z [V/m] @ f=',num2str(f)]);caxis(([-max(max(abs(PEz))) max(max(abs(PEz)))]));colorbar;
    figure(); pcolor(zz_mesh,rr_mesh,(abs(PEphi))); shading interp; hold on; contour(zz_mesh,rr_mesh,(abs(PEphi))); hold off;
    xlabel('Length [m]'); ylabel('Radius [m]'); title(['Field E_{\phi} [V/m] @ f=',num2str(f)]);caxis([-max(max(abs(PEphi))) max(max(abs(PEphi)))]);colorbar;
    figure(); pcolor(zz_mesh,rr_mesh,(abs(PHr))); shading interp; hold on; contour(zz_mesh,rr_mesh,(abs(PHr))); hold off;
    xlabel('Length [m]'); ylabel('Radius [m]'); title(['Field H_r [A/m] @ f=',num2str(f)]); caxis([-max(max(abs(PHr))) max(max(abs(PHr)))]); colorbar;
    figure(); pcolor(zz_mesh,rr_mesh,(abs(PHz))); shading interp; hold on; contour(zz_mesh,rr_mesh,(abs(PHz))); hold off;
    xlabel('Length [m]'); ylabel('Radius [m]'); title(['Field H_z [A/m] @ f=',num2str(f)]);caxis([-max(max(abs(PHz))) max(max(abs(PHz)))]);colorbar;
    figure(); pcolor(zz_mesh,rr_mesh,(abs(PHphi))); shading interp; hold on; contour(zz_mesh,rr_mesh,(abs(PHphi))); hold off;
    xlabel('Length [m]'); ylabel('Radius [m]'); title(['Field H_{\phi} [A/m] @ f=',num2str(f)]);caxis([-max(max(abs(PHphi))) max(max(abs(PHphi)))]);colorbar;
   

%% Write fields
SaveDir='~/scratch0/FiniteLength/ModeMatching/Transverse_simulations/';  
material={};
material.name='test';
    nome_fileEz=['MMMdip_Ez_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_phi',num2str(ph),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];
    nome_fileEr=['MMMdip_Er_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_phi',num2str(ph),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];
    nome_fileEphi=['MMMdip_Ephi_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_phi',num2str(ph),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];
    nome_fileHz=['MMMdip_Hz_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_phi',num2str(ph),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];
    nome_fileHr=['MMMdip_Hr_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_phi',num2str(ph),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];
    nome_fileHphi=['MMMdip_Hphi_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_phi',num2str(ph),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];
    nome_fileMeshGridZ=['MMMdip_MeshGridZ_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_phi',num2str(ph),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];
    nome_fileMeshGridR=['MMMdip_MeshGridR_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_phi',num2str(ph),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];   
    dlmwrite([SaveDir,nome_fileEz],PEz,'precision',10);
    dlmwrite([SaveDir,nome_fileEr],PEr,'precision',10);
    dlmwrite([SaveDir,nome_fileEphi],PEphi,'precision',10);
    dlmwrite([SaveDir,nome_fileHz],PHz,'precision',10);
    dlmwrite([SaveDir,nome_fileHr],PHr,'precision',10);
    dlmwrite([SaveDir,nome_fileHphi],PHphi,'precision',10);
    dlmwrite([SaveDir,nome_fileMeshGridZ],[zz_mesh],'precision',10);
    dlmwrite([SaveDir,nome_fileMeshGridR],[rr_mesh],'precision',10);
    disp([nome_fileEz]);
    disp([nome_fileEr]);
    disp([nome_fileEphi]);
    disp([nome_fileHz]);
    disp([nome_fileHr]);
    disp([nome_fileHphi]);
    disp([nome_fileMeshGridR]);
    disp([nome_fileMeshGridZ]);
    end
end

%% Save
%  Define files
nome_fileRe=['MMMdip_Re_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_Material_',material.name,...
         '_fmin',num2str(min(f_vec)),...
         '_fmax',num2str(max(f_vec)),'_P',num2str(P),'_S',num2str(S),'.txt'];
nome_fileIm=['MMMdip_Im_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_Material_',material.name,...
         '_fmin',num2str(min(f_vec)),...
         '_fmax',num2str(max(f_vec)),'_P',num2str(P),'_S',num2str(S),'.txt'];
nome_filePP=['MMMdip_Postproc_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t','_Material_',material.name,...
         '_fmin',num2str(min(f_vec)),...
         '_fmax',num2str(max(f_vec)),'_P',num2str(P),'_S',num2str(S),'.txt'];  
% Write impedance
dlmwrite([SaveDir,nome_fileRe],[f_vec',real(Z)'],'precision',10);
dlmwrite([SaveDir,nome_fileIm],[f_vec',imag(Z)'],'precision',10);
disp(['Data saved in:',SaveDir]);
disp([nome_fileRe]);
disp([nome_fileIm]);

end