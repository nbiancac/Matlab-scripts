function [freq, Z,TMmodes_WG,TMmodes_C]=MMM_longitudinal(beta,b,t,L,material,fin,fstep,fout,fdiscrete,f_post, P,S,SaveDir,fields,loss)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               MMM_longitudinal.m
%   
%   Impedance of a Finite length device with the Mode Matching Method (v2.05.2014). 
%
%   function [freq, Z,TMmodes_WG,TMmodes_C] = MMM_longitudinal(beta,b,t,L,material,fin,fstep,fout,fdiscrete,f_post, P,S,SaveDir,fields,loss)
%   
%   Function description:
%     This function calculates the longitudinal coupling impedance for a
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
%     fields    :       Calculates fields pattern at each frequency in f_post;
%     loss      :       Calculate power loss at each frequency in  f_post;
%
%   Output:
%     freq      :       Frequencies used [Hz]
%     Zlong     :       Longitudinal impedance [Ohm] 
%     TMmodes_WG:       TM modes frequencies in the pipes [GHz]  
%     TMmodes_C :       TM modes frequencies in the empty cavity [GHz]
%  
%   Written Files:
%     "MMMlong_Re_*.txt": Real part of Zlong;
%     "MMMlong_Im_*.txt": Imaginary part of Zlong;
%     "MMMlong_Postproc_*.txt": File with postprocessing results. 
%                       It contains the postprocessed frequencies and the power losses for radiation,
%                       volume conductivity, volume electric and magnetic losses.
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
% [freq,Zlong,~,~]=MMM_longitudinal(beta,b,t,L,material,fin,fstep,fout,fdiscrete,f_post,P,S,SaveDir,0,0);
%     
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Physical constants
sl         =       299792448.2;      % speed of light [m/s]
epsilon0   =       10^7/(4*pi*sl^2); % Vacuum dielectric constant [Farad/m]
mu0        =       4*pi*10^-7;       % Vacuum permeability [Henry^-1/m]
Z0         =       120*pi;           % Vacuum characteristic impedance [Ohm]
v          =       beta*sl;
%% Beam
Q       =       1;                % Beam charge [C]
gamma   =       1/sqrt(1-beta^2); % Relativistic gamma

%% Number of modes
    
nu=0;
% [Jzero,Jzeroprime]=besselzeros2(nu+1,P);
% Jzeroprime=Jzeroprime(:,nu+1);
% Jzero=Jzero(:,nu+1); % Zeros in bessel J0(x) 
Jzero=zerobess('J',nu,P);
Jzeroprime=zerobess('DJ',nu,P);

%% Matrix inizialization
    
A_TM_vec=[];
C_TM_vec=[];
D_TM_vec=[];
W_TM_smaterial.freq=[];
H_TM_smaterial.freq=[];
Z_sx_vec=[];
Z_dx_vec=[];
Z_cav_vec=[];
Zrew=[];
material.freq=[];
B_vec=[];


Z=[];
B =zeros(S,1);
IIss=zeros(S,S);
II=zeros(S,S);
M1=zeros(S,P);
M2=zeros(P,S);
N1=zeros(S,S);
N2=zeros(P,P);
U1=zeros(S,S);
Z_integral_sx=zeros(1,P);
Z_integral_dx=zeros(1,P);
Z_C=zeros(1,P);
Z_D=zeros(1,P);
Z_A=zeros(1,S);
Z=zeros(1,length(material.freq));

A_TM_post=[];
C_TM_post=[];
D_TM_post=[];    
%% Dimensions

c = b+t;     % cavity radius

%% Frequency scan

if ~isempty(fdiscrete)
    material.freq=fdiscrete;
    if size(fdiscrete,1)>size(fdiscrete,2) material.freq=material.freq'; end
else
    material.freq=fin:fstep:fout; 
end
if size(f_post,1)>size(f_post,2) f_post=f_post'; end

material.freq=unique(sort([material.freq,f_post]));
freq=material.freq;

for ff=1:length(material.freq)
    
    
    f=freq(ff);
    disp(f);
    
    %% Material properties
    
    omega=2*pi*f;
    if ischar(material.e_r)
        epsilonf=epsilon0*(eval(material.e_r)-1i*material.sigma/(omega*epsilon0)); % Relative permeability
    elseif length(material.e_r)==length(material.freq)
        epsilonf=epsilon0*((material.e_r(ff))-1i*material.sigma/(omega*epsilon0)); 
    end
    if ischar(material.mu_r)    % Relative susceptibility
        muf=mu0*eval(material.mu_r);
    elseif length(material.mu_r)==length(material.freq)
        muf=mu0*material.mu_r(ff);
    end
    
    Zf=sqrt(muf/epsilonf);
    Yf=1/Zf;
    
    
    %% Kappas
    k0=omega*sqrt(mu0*epsilon0);
    kf=omega*sqrt(muf*epsilonf);
    kb=omega/(beta*sl);

    for s=0:S-1;
        
        ss=s+1; %indici per riempire le matrici
        if s==0 es=1;  else es=2;    end
        if s==0 delta_s0=1; else delta_s0=0; end
           
        
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
            alpha_ps2         =       alpha_p^2+alpha_s^2;
            alpha_ps          =       sqrt(alpha_ps2);
            alpha_tilde_s2    =       alpha_0^2-alpha_s^2;
            alpha_tilde_p2    =       alpha_0^2-alpha_p^2;
            alpha_tilde_sf2   =       alpha_f^2-alpha_s^2;
            alpha_tilde_p     =       conj(sqrt(alpha_tilde_p2));
            alpha_tilde_s     =       (sqrt(alpha_tilde_s2));
            alpha_tilde_sf    =       sqrt(alpha_tilde_sf2);
            
            if material.sigma==0 && imag(epsilonf)==0
                
            alpha_tilde_sf               =        conj(alpha_tilde_sf);
            end
            
            
            in =        alpha_tilde_sf;
           
            W_TM_sf            =        1-(besselh(0,2,in*c/b,1))/(besselh(0,1,in*c/b,1))*besselh(0,1,in,1)/besselh(0,2,in,1)*exp(-2*1i*in*t/b);
            W_TMprime_sf       =        -besselh(1,2,in,1)/besselh(0,2,in,1)+(besselh(1,1,in,1))*(besselh(0,2,in*c/b,1))/besselh(0,1,in*c/b,1)/besselh(0,2,in,1)*exp(-2*1i*in*t/b);
            H_TM_sf            =        (1i*alpha_f*sqrt(L/es)*W_TMprime_sf)/(Zf*alpha_tilde_sf) + (1i*besselj(1, alpha_tilde_s,1)*alpha_0*sqrt(L/es)*W_TM_sf)/(besselj(0, alpha_tilde_s,1)*Z0*alpha_tilde_s); 
            
            % Matrices for coefficients
            
            B(ss,1)            =        ((-1i/2)*(exp(1i*pi*s) - exp((1i*L*alpha_b)/b))*alpha_b)/(exp((1i*L*alpha_b)/b)*pi*besseli(0, alpha_b/gamma)*(alpha_b^2 - alpha_s^2)*H_TM_sf);
            IIss(ss,ss)        =        (-1)^s;
            II(ss,ss)          =        1;
            M1(ss,pp)          =        (sign(besselj(1, alpha_p))*alpha_tilde_p)/(alpha_p^2*(alpha_02 - alpha_ps^2));
            M2(pp,ss)          =        ((exp(((-1i)*L*alpha_tilde_p)/b)*(-1)^s - 1)*L*W_TM_sf)/((alpha_02 - alpha_ps^2)*sqrt(L/es));
            N1(ss,ss)          =        (b*alpha_0)/(sqrt(pi)*Z0*H_TM_sf);
            N2(pp,pp)          =        (sqrt(pi)*sign(besselj(1, alpha_p))*alpha_p^2)/L;
                
            % Matrices for impedance
            Z_integral_sx(1,pp)   =    (1i*b)/(sqrt(pi)*abs(besselj(1, alpha_p))*alpha_p*(alpha_b + alpha_tilde_p));
            Z_integral_dx(1,pp)   =    ((-1i)*b*exp((1i*L*alpha_b)/b))/(sqrt(pi)*abs(besselj(1, alpha_p))*alpha_p*(alpha_b - alpha_tilde_p));
%             Z_C(1,pp)             =    (b*((-exp((1i*L*alpha_b)/b) + cos((L*alpha_tilde_p)/b))*csc((L*alpha_tilde_p)/b)*alpha_b + 1i*alpha_tilde_p))/(sqrt(pi)*abs(besselj(1, alpha_p))*alpha_p*(alpha_b^2 - alpha_tilde_p2));
%             Z_D(1,pp)             =    -((b*((-1 + exp((1i*L*alpha_b)/b)*cos((L*alpha_tilde_p)/b))*csc((L*alpha_tilde_p)/b)*alpha_b - 1i*exp((1i*L*alpha_b)/b)*alpha_tilde_p))/(sqrt(pi)*abs(besselj(1, alpha_p))*alpha_p*(alpha_b^2 - alpha_tilde_p2)));
            Z_C(1,pp)             =    (b*((-2i*exp(1i*(L*alpha_b/b-L*alpha_tilde_p/b))/(1-exp(-2i*L*alpha_tilde_p/b)) + cot((L*alpha_tilde_p)/b))*alpha_b + 1i*alpha_tilde_p))/(sqrt(pi)*abs(besselj(1, alpha_p))*alpha_p*(alpha_b^2 - alpha_tilde_p2));
            Z_D(1,pp)             =    -((b*((-2i*exp(-1i*L*alpha_tilde_p/b)/(1-exp(-2i*L*alpha_tilde_p/b)) + exp((1i*L*alpha_b)/b)*cot((L*alpha_tilde_p)/b))*alpha_b - 1i*exp((1i*L*alpha_b)/b)*alpha_tilde_p))/(sqrt(pi)*abs(besselj(1, alpha_p))*alpha_p*(alpha_b^2 - alpha_tilde_p2)));
            Z_A(1,ss)             =    (1i*b*(-1 + (-1)^s*exp((1i*L*alpha_b)/b))*alpha_b*W_TM_sf)/(besselj(0, alpha_tilde_s)*(alpha_b^2 - alpha_s^2)*sqrt(L/es));
        end
        
    end
    
          
    
    % Inversions
  
    
    U1=II - N1*M1*N2*M2 - N1*IIss*M1*N2*M2*IIss;
   
    A_TM=U1\B;
    D_TM=N2*M2*IIss*A_TM;
    C_TM=N2*M2*A_TM;
    
    Zsx=(Z_integral_sx)*C_TM;
    Zdx=(Z_integral_dx)*D_TM;
    Zcav=Z_C*C_TM+Z_D*D_TM+Z_A*A_TM;
    Z(ff)=Zsx+Zdx+Zcav;
    A_TM_vec=[A_TM_vec,A_TM];
    C_TM_vec=[C_TM_vec,C_TM];
    D_TM_vec=[D_TM_vec,D_TM];
    W_TM_smaterial.freq=[W_TM_smaterial.freq,W_TM_sf];
    H_TM_smaterial.freq=[H_TM_smaterial.freq,H_TM_sf];
    B_vec=[B_vec,B];
    
    Z_sx_vec=[Z_sx_vec,Zsx];
    Z_dx_vec=[Z_dx_vec,Zdx];
    Z_cav_vec=[Z_cav_vec,Zcav];
  
    if ismember(f,f_post)
        A_TM_post=[A_TM_post,A_TM];
        C_TM_post=[C_TM_post,C_TM];
        D_TM_post=[D_TM_post,D_TM];
        
    end

end

TMmodes_WG=Jzero*sl/2/pi/b/1e9;
% disp(TMmodes_WG);
TMmodes_C=[];
for s=0:S-1
    TMmodes_C=[TMmodes_C,sl/2/pi/c/1e9*sqrt(Jzero.^2+(s*pi/L*c).^2)];
end
% disp(TMmodes_C)
%% Plottery

% figure(1);
%     subplot(211) 
%     plot(material.freq,real(Z),'b'); hold on;
%     grid on;
%     title('Real part')
%     ylabel('Zlong [Ohm]','FontSize',12);
%     xlabel('frequency (Hz)','FontSize',12);
%     legend('Mode Matching')
%        
%     subplot(212) 
%     plot(material.freq,imag(Z),'r');hold on;
%     grid on;
%     title('Imaginary part')
%     xlabel('frequency (Hz)','FontSize',12);
%     ylabel('Zlong [Ohm]','FontSize',12);
%     legend('Mode Matching')

%% Postprocessing
% Power dissipated
W1=[];W2=[];W3=[];W4=[];W5=[];W6=[];W7=[];Qvec=[];
Fields_vec=[];


for f=f_post
    if loss==1                     
    disp(['Postprocessing power dissipated at ',num2str(f),'...']);
    
    % Material properties
    
    omega=2*pi*f;
    epsilonf=epsilon0*(eval(material.e_r)-1i*material.sigma/(omega*epsilon0)); % Relative dielectric 
    muf=mu0*eval(material.mu_r);
    Zf=sqrt(muf/epsilonf);
    Yf=1/Zf;
    
    % Kappas
    
    k0=omega*sqrt(mu0*epsilon0);
    kf=omega*sqrt(muf*epsilonf);
    kb=omega/(beta*sl);
  
    % Meshing
    
    Dzz=L/15;
    Drr=t/15;
    P1=[];P2=[];P3=[];P4=[];P5=[];P6=[];P7=[];
    PEz=[];PEr=[];PHphi=[]; PPoyr=[];PPoyz=[];
    rr_mesh=b:Drr:c;
    zz_mesh=0:Dzz:L;
    
    % Coefficients
    
    [~,indf]=ismember(f,f_post);
    A_TM=A_TM_post(:,indf);
    
    for r_ind=1:length(rr_mesh)
        rr=rr_mesh(r_ind);
        
       for z_ind=1:length(zz_mesh)
           zz=zz_mesh(z_ind);
            
            
                
                Ez=0;
                Er=0;
                Hphi=0;
               
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
                    alpha_tilde_s     =       sqrt(alpha_tilde_s2);
                    alpha_tilde_sf    =       sqrt(alpha_tilde_sf2);

                    if material.sigma==0
                        alpha_tilde_sf    =    conj(alpha_tilde_sf);
                    end

                    in                 =        alpha_tilde_sf;

                    W_TM_sf            =        1-(besselh(0,2,in*c/b,1))/(besselh(0,1,in*c/b,1))*besselh(0,1,in*rr/b,1)/besselh(0,2,in*rr/b,1)*exp(-2*1i*in*(c-rr)/b);
                    W_TMprime_sf       =        -besselh(1,2,in*rr/b,1)/besselh(0,2,in*rr/b,1)+(besselh(1,1,in*rr/b,1))*(besselh(0,2,in*c/b,1))/besselh(0,1,in*c/b,1)/besselh(0,2,in*rr/b,1)*exp(-2*1i*in*(c-rr)/b);
                   
                    % field coefficients
                   Ez=Ez+cos(zz*alpha_s/b)*A_TM(ss)*besselh(0,2,in*rr/b)*W_TM_sf/(sqrt(L/es));  % Ez
                   Er=Er-sin(zz*alpha_s/b)*A_TM(ss)*besselh(0,2,in*rr/b)*alpha_s*W_TMprime_sf/(sqrt(L/es)*alpha_tilde_sf); % Er
                   Hphi=Hphi+((-1i)*A_TM(ss)*besselh(0,2,in*rr/b)*cos((zz*alpha_s)/b)*alpha_f*W_TMprime_sf)/(Zf*sqrt(L/es)*alpha_tilde_sf); % Hphi
         
                    
                end
                
            if rr==b
            P1(z_ind)=(0.5)*b*Ez*conj(Hphi);
            end
            P2(r_ind,z_ind)=1/2*material.sigma*(Ez*conj(Ez)+Er*conj(Er));
            P3(r_ind,z_ind)=1/2*1i*(mu0*eval(material.mu_r))*omega*Hphi*conj(Hphi);
            P4(r_ind,z_ind)=-1/2*1i*(epsilon0*eval(material.e_r))*omega*(Er*conj(Er)+Ez*conj(Ez));
            PEz(r_ind,z_ind)=Ez;
            PEr(r_ind,z_ind)=Er;
            PHphi(r_ind,z_ind)=Hphi;
            PPoyr(r_ind,z_ind)=-0.5*Ez*conj(Hphi);
            PPoyz(r_ind,z_ind)=0.5*Er*conj(Hphi);
            
      end % z loop
  
    end % r loop
    
    W1=[W1,2*pi*trapz(zz_mesh,P1,2)];
    W2=[W2,2*pi*trapz(rr_mesh,rr_mesh.*conj(trapz(zz_mesh,P2,2)'))];
    W3=[W3,2*pi*trapz(rr_mesh,rr_mesh.*conj(trapz(zz_mesh,P3,2)'))];
    W4=[W4,2*pi*trapz(rr_mesh,rr_mesh.*conj(trapz(zz_mesh,P4,2)'))];
    
    
    % Meshing for cavity losses
    
    Dzz=L/15;
    Drr=b/15;
   
    PEz=[];
    PEr=[];
    PHphi=[];
    rr_mesh=[Drr/2:Drr:c];
    zz_mesh=[-L-Dzz/2:Dzz:2*L];
    
    % Coefficients
    
    [~,indf]=ismember(f,f_post);
    A_TM=A_TM_post(:,indf);
    C_TM=C_TM_post(:,indf);
    D_TM=D_TM_post(:,indf);
    
    for r_ind=1:length(rr_mesh)
        rr=rr_mesh(r_ind);
        
       for z_ind=1:length(zz_mesh)
           zz=zz_mesh(z_ind);
                 
                Ez=0;
                Er=0;
                Hphi=0;
                Ez_b=0;
                Er_b=0;
                Hphi_b=0;
                % CAVITY
                P5(r_ind,z_ind)=0;
                P6(r_ind,z_ind)=0;
                P7(r_ind,z_ind)=0;
                
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
                            alpha_ps2         =       alpha_p^2+alpha_s^2;
                            alpha_ps          =       sqrt(alpha_ps2);
                            alpha_tilde_s2    =       alpha_0^2-alpha_s^2;
                            alpha_tilde_p2    =       alpha_0^2-alpha_p^2;
                            alpha_tilde_sf2   =       alpha_f^2-alpha_s^2;
                            alpha_tilde_p     =       conj(sqrt(alpha_tilde_p2));
                            alpha_tilde_s     =       (sqrt(alpha_tilde_s2));
                            alpha_tilde_sf    =       sqrt(alpha_tilde_sf2);

                            if material.sigma==0
                            alpha_tilde_sf               =        conj(alpha_tilde_sf);
                            end

                            in                 =        alpha_tilde_sf;

                            W_TM_sf            =        1-(besselh(0,2,in*c/b,1))/(besselh(0,1,in*c/b,1))*besselh(0,1,in,1)/besselh(0,2,in,1)*exp(-2*1i*in*t/b);
                            W_TMprime_sf       =        -besselh(1,2,in,1)/besselh(0,2,in,1)+(besselh(1,1,in,1))*(besselh(0,2,in*c/b,1))/besselh(0,1,in*c/b,1)/besselh(0,2,in,1)*exp(-2*1i*in*t/b);
                            H_TM_sf            =        (1i*alpha_f*sqrt(L/es)*W_TMprime_sf)/(Zf*alpha_tilde_sf) + (1i*besselj(1, alpha_tilde_s)*alpha_0*sqrt(L/es)*W_TM_sf)/(besselj(0, alpha_tilde_s)*Z0*alpha_tilde_s); 

                            I_TM=(1i*b*alpha_0*((1i*b*(C_TM(pp))*alpha_tilde_p)/(alpha_p^2*sqrt(L/es)) + (1i*(-1)^s*b*(D_TM(pp))*alpha_tilde_p)/(alpha_p^2*sqrt(L/es)) - 2*sqrt(pi)*sign(besselj(1, alpha_p))*A_TM(ss)*besselh(0,2,in)*W_TM_sf))/(Z0*(alpha_0^2 - alpha_ps^2));
                        
                            Ez=Ez+I_TM*(((-1i)*Z0*alpha_ps)/alpha_0)*(besselj(0, (rr*alpha_p)/b)*cos((zz*alpha_s)/b)*alpha_p)/(b*sqrt(pi)*abs(besselj(1, alpha_p))*alpha_ps*sqrt(L/es));
                            Er=Er+-I_TM*(((-1i)*Z0*alpha_ps)/alpha_0)*((besselj(1, (rr*alpha_p)/b)*sin((zz*alpha_s)/b)*alpha_s)/(b*sqrt(pi)*abs(besselj(1, alpha_p))*alpha_ps*sqrt(L/es)));
                            Hphi=Hphi+I_TM*(besselj(1, (rr*alpha_p)/b)*cos((zz*alpha_s)/b))/(b*sqrt(pi)*abs(besselj(1, alpha_p))*sqrt(L/es));
                        
                       end
                    end
                    
                
                end
                
                P5(r_ind,z_ind)=1/2*0*(Ez*conj(Ez)+Er*conj(Er));
                P6(r_ind,z_ind)=1/2*1i*(mu0)*omega*Hphi*conj(Hphi);
                P7(r_ind,z_ind)=-1/2*1i*(epsilon0)*omega*(Er*conj(Er)+Ez*conj(Ez));
                                
       end % end z loop
       
            
    end % end r loop
    
    
    W5=[W5,2*pi*trapz(rr_mesh,rr_mesh.*conj(trapz(zz_mesh,P5,2)'))];
    W6=[W6,2*pi*trapz(rr_mesh,rr_mesh.*conj(trapz(zz_mesh,P6,2)'))];
    W7=[W7,2*pi*trapz(rr_mesh,rr_mesh.*conj(trapz(zz_mesh,P7,2)'))];
    Qvec=[Qvec,imag(W6(end)+W3(end))/W2(end)];
    end

    if fields==1
    disp(['Postprocessing Fields  at ',num2str(f),'...']);
        
    % Material properties
    
    omega=2*pi*f;
    epsf=epsilon0*eval(material.e_r);
    epsilonf=epsilon0*(eval(material.e_r)-1i*material.sigma/(omega*epsilon0)); % Relative dielectric 
    muf=mu0*eval(material.mu_r);
    Zf=sqrt(muf/epsilonf);
    Yf=1/Zf;
    
    % Kappas
    
    k0=omega*sqrt(mu0*epsilon0);
    kf=omega*sqrt(muf*epsilonf);
    kb=omega/(beta*sl);
  
    % Meshing
    
    Dzz=L/15;
    Drr=b/15;
   
    PEz=[];
    PEr=[];
    PHphi=[];
    rr_mesh=[Drr/2:Drr:c];
    zz_mesh=[-L-Dzz/2:Dzz:2*L];
    
    % Coefficients
    
    [~,indf]=ismember(f,f_post);
    A_TM=A_TM_post(:,indf);
    C_TM=C_TM_post(:,indf);
    D_TM=D_TM_post(:,indf);
    
    for r_ind=1:length(rr_mesh)
        rr=rr_mesh(r_ind);
        
       for z_ind=1:length(zz_mesh)
           zz=zz_mesh(z_ind);
                 
                Ez=0;
                Er=0;
                Hphi=0;
                Ez_b=0;
                Er_b=0;
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
                        alpha_tilde_s     =       sqrt(alpha_tilde_s2);
                        alpha_tilde_sf    =       sqrt(alpha_tilde_sf2);

                        if material.sigma==0
                            alpha_tilde_sf    =    conj(alpha_tilde_sf);
                        end

                        in                 =        alpha_tilde_sf;

                        W_TM_sf            =        1-(besselh(0,2,in*c/b,1))/(besselh(0,1,in*c/b,1))*besselh(0,1,in*rr/b,1)/besselh(0,2,in*rr/b,1)*exp(-2*1i*in*(c-rr)/b);
                        W_TMprime_sf       =        -besselh(1,2,in*rr/b,1)/besselh(0,2,in*rr/b,1)+(besselh(1,1,in*rr/b,1))*(besselh(0,2,in*c/b,1))/besselh(0,1,in*c/b,1)/besselh(0,2,in*rr/b,1)*exp(-2*1i*in*(c-rr)/b);

                        Ez=Ez     +A_TM(ss)*cos(zz*alpha_s/b)*besselh(0,2,in*rr/b)*W_TM_sf/(sqrt(L/es));  % Ez
                        Er=Er     -A_TM(ss)*sin(zz*alpha_s/b)*besselh(0,2,in*rr/b)*alpha_s*W_TMprime_sf/(sqrt(L/es)*alpha_tilde_sf); % Er
                        Hphi=Hphi +A_TM(ss)*((-1i)*besselh(0,2,in*rr/b)*cos((zz*alpha_s)/b)*alpha_f*W_TMprime_sf)/(Zf*sqrt(L/es)*alpha_tilde_sf); % Hphi
                    
                        Ez_b=0;
                        Er_b=0;
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
                        alpha_f           =       kf*b;
                        alpha_f2          =       alpha_f^2;
                        alpha_b           =       kb*b;
                        alpha_b2          =       alpha_b^2;
                        alpha_p           =       Jzero(p);
                        alpha_p2          =       alpha_p^2;
                        alpha_ps2         =       alpha_p^2+alpha_s^2;
                        alpha_ps          =       sqrt(alpha_ps2);
                        alpha_tilde_p2    =       alpha_0^2-alpha_p^2;
                        alpha_tilde_p     =       conj(sqrt(alpha_tilde_p2));

                        Ez=Ez    +C_TM(pp)*(exp((1i*zz*alpha_tilde_p)/b)*besselj(0, (rr*alpha_p)/b))/(sqrt(pi)*abs(besselj(1, alpha_p))*alpha_p);  % Ez
                        Er=Er    +C_TM(pp)*((-1i)*exp((1i*zz*alpha_tilde_p)/b)*besselj(1, (rr*alpha_p)/b)*alpha_tilde_p)/(sqrt(pi)*abs(besselj(1, alpha_p))*alpha_p^2); % Er
                        Hphi=Hphi+C_TM(pp)*(1i*exp((1i*zz*alpha_tilde_p)/b)*besselj(1, (rr*alpha_p)/b)*alpha_0)/(sqrt(pi)*abs(besselj(1, alpha_p))*Z0*alpha_p^2); % Hphi
                   
                        Er_b=-(omega*(-((besseli(1, (rr*alpha_b)/(b*gamma))*besselk(0, alpha_b/gamma))/besseli(0, alpha_b/gamma)) - besselk(1, (rr*alpha_b)/(b*gamma))))/(2*exp((1i*zz*alpha_b)/b)*pi*v^2*gamma*epsilon0);
                        Ez_b=((1i/2)*omega*(-((besseli(0, (rr*alpha_b)/(b*gamma))*besselk(0, alpha_b/gamma))/besseli(0, alpha_b/gamma)) + besselk(0, (rr*alpha_b)/(b*gamma))))/(exp((1i*zz*alpha_b)/b)*pi*v^2*gamma^2*epsilon0);
                        Hphi_b=-(omega*(-((besseli(1, (rr*alpha_b)/(b*gamma))*besselk(0, alpha_b/gamma))/besseli(0, alpha_b/gamma)) - besselk(1, (rr*alpha_b)/(b*gamma))))/(2*exp((1i*zz*alpha_b)/b)*pi*v*gamma);
                        
                    end
               
                end
                
                % PIPE DX
                
                if rr<b && zz>L
               
                    for p=1:P
                        pp=p;
                        alpha_0           =       k0*b;
                        alpha_02          =       alpha_0^2;
                        alpha_f           =       kf*b;
                        alpha_f2          =       alpha_f^2;
                        alpha_b           =       kb*b;
                        alpha_b2          =       alpha_b^2;
                        alpha_p           =       Jzero(p);
                        alpha_p2          =       alpha_p^2;
                        alpha_ps2         =       alpha_p^2+alpha_s^2;
                        alpha_ps          =       sqrt(alpha_ps2);
                        alpha_tilde_p2    =       alpha_0^2-alpha_p^2;
                        alpha_tilde_p     =       conj(sqrt(alpha_tilde_p2));

                        Ez=Ez+D_TM(pp)*besselj(0, (rr*alpha_p)/b)/(exp((1i*(-L + zz)*alpha_tilde_p)/b)*sqrt(pi)*abs(besselj(1, alpha_p))*alpha_p);  % Ez
                        Er=Er+D_TM(pp)*(1i*besselj(1, (rr*alpha_p)/b)*alpha_tilde_p)/(exp((1i*(-L + zz)*alpha_tilde_p)/b)*sqrt(pi)*abs(besselj(1, alpha_p))*alpha_p^2); % Er
                        Hphi=Hphi+D_TM(pp)*(1i*besselj(1, (rr*alpha_p)/b)*alpha_0)/(exp((1i*(-L + zz)*alpha_tilde_p)/b)*sqrt(pi)*abs(besselj(1, alpha_p))*Z0*alpha_p^2); % Hphi
                        
                        Er_b=-(omega*(-((besseli(1, (rr*alpha_b)/(b*gamma))*besselk(0, alpha_b/gamma))/besseli(0, alpha_b/gamma)) - besselk(1, (rr*alpha_b)/(b*gamma))))/(2*exp((1i*zz*alpha_b)/b)*pi*v^2*gamma*epsilon0);
                        Ez_b=((1i/2)*omega*(-((besseli(0, (rr*alpha_b)/(b*gamma))*besselk(0, alpha_b/gamma))/besseli(0, alpha_b/gamma)) + besselk(0, (rr*alpha_b)/(b*gamma))))/(exp((1i*zz*alpha_b)/b)*pi*v^2*gamma^2*epsilon0);
                        Hphi_b=-(omega*(-((besseli(1, (rr*alpha_b)/(b*gamma))*besselk(0, alpha_b/gamma))/besseli(0, alpha_b/gamma)) - besselk(1, (rr*alpha_b)/(b*gamma))))/(2*exp((1i*zz*alpha_b)/b)*pi*v*gamma);
                        
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
                            alpha_ps2         =       alpha_p^2+alpha_s^2;
                            alpha_ps          =       sqrt(alpha_ps2);
                            alpha_tilde_s2    =       alpha_0^2-alpha_s^2;
                            alpha_tilde_p2    =       alpha_0^2-alpha_p^2;
                            alpha_tilde_sf2   =       alpha_f^2-alpha_s^2;
                            alpha_tilde_p     =       conj(sqrt(alpha_tilde_p2));
                            alpha_tilde_s     =       (sqrt(alpha_tilde_s2));
                            alpha_tilde_sf    =       sqrt(alpha_tilde_sf2);

                            if material.sigma==0
                            alpha_tilde_sf               =        conj(alpha_tilde_sf);
                            end

                            in                 =        alpha_tilde_sf;

                            W_TM_sf            =        1-(besselh(0,2,in*c/b,1))/(besselh(0,1,in*c/b,1))*besselh(0,1,in,1)/besselh(0,2,in,1)*exp(-2*1i*in*t/b);
                            W_TMprime_sf       =        -besselh(1,2,in,1)/besselh(0,2,in,1)+(besselh(1,1,in,1))*(besselh(0,2,in*c/b,1))/besselh(0,1,in*c/b,1)/besselh(0,2,in,1)*exp(-2*1i*in*t/b);
                            H_TM_sf            =        (1i*alpha_f*sqrt(L/es)*W_TMprime_sf)/(Zf*alpha_tilde_sf) + (1i*besselj(1, alpha_tilde_s)*alpha_0*sqrt(L/es)*W_TM_sf)/(besselj(0, alpha_tilde_s)*Z0*alpha_tilde_s); 

                            I_TM=(1i*b*alpha_0*((1i*b*(C_TM(pp))*alpha_tilde_p)/(alpha_p^2*sqrt(L/es)) + (1i*(-1)^s*b*(D_TM(pp))*alpha_tilde_p)/(alpha_p^2*sqrt(L/es)) - 2*sqrt(pi)*sign(besselj(1, alpha_p))*A_TM(ss)*besselh(0,2,in)*W_TM_sf))/(Z0*(alpha_0^2 - alpha_ps^2));
                        
                            Ez=Ez+I_TM*(((-1i)*Z0*alpha_ps)/alpha_0)*(besselj(0, (rr*alpha_p)/b)*cos((zz*alpha_s)/b)*alpha_p)/(b*sqrt(pi)*abs(besselj(1, alpha_p))*alpha_ps*sqrt(L/es));
                            Er=Er+-I_TM*(((-1i)*Z0*alpha_ps)/alpha_0)*((besselj(1, (rr*alpha_p)/b)*sin((zz*alpha_s)/b)*alpha_s)/(b*sqrt(pi)*abs(besselj(1, alpha_p))*alpha_ps*sqrt(L/es)));
                            Hphi=Hphi+I_TM*(besselj(1, (rr*alpha_p)/b)*cos((zz*alpha_s)/b))/(b*sqrt(pi)*abs(besselj(1, alpha_p))*sqrt(L/es));
                        
                            Er_b=-(omega*(-((besseli(1, (rr*alpha_b)/(b*gamma))*besselk(0, alpha_b/gamma))/besseli(0, alpha_b/gamma)) - besselk(1, (rr*alpha_b)/(b*gamma))))/(2*exp((1i*zz*alpha_b)/b)*pi*v^2*gamma*epsilon0);
                            Ez_b=((1i/2)*omega*(-((besseli(0, (rr*alpha_b)/(b*gamma))*besselk(0, alpha_b/gamma))/besseli(0, alpha_b/gamma)) + besselk(0, (rr*alpha_b)/(b*gamma))))/(exp((1i*zz*alpha_b)/b)*pi*v^2*gamma^2*epsilon0);
                            Hphi_b=-(omega*(-((besseli(1, (rr*alpha_b)/(b*gamma))*besselk(0, alpha_b/gamma))/besseli(0, alpha_b/gamma)) - besselk(1, (rr*alpha_b)/(b*gamma))))/(2*exp((1i*zz*alpha_b)/b)*pi*v*gamma);
                        end
                    end
                    
                    
                    
                
                end
                
                if rr>b
                    if zz>L || zz<0
                           Ez=nan;
                           Er=nan;
                           Hphi=nan;
                           Ez_b=0;
                           Er_b=0;
                           Hphi_b=0;
                    end
                end
                
                PEz(r_ind,z_ind)=Ez+Ez_b; % Field vectors
                PEr(r_ind,z_ind)=Er+Er_b;
                PHphi(r_ind,z_ind)=Hphi+1*Hphi_b; 
       
       end % end z loop
       
            
    end % end r loop
      
    
    % Plot Ez,r Hphi
   
%     figure(); pcolor(zz_mesh,rr_mesh,(real(PEr))); shading interp; hold on; contour(zz_mesh,rr_mesh,(real(PEr))); hold off;
%     xlabel('Length [m]'); ylabel('Radius [m]'); title('Field E_r [V/m]'); caxis([-max(max(real(PEr))) max(max(real(PEr)))]); colorbar;
%     figure(); pcolor(zz_mesh,rr_mesh,(real(PEz))); shading interp; hold on; contour(zz_mesh,rr_mesh,(real(PEz))); hold off;
%     xlabel('Length [m]'); ylabel('Radius [m]'); title('Field E_z [V/m]');caxis([-max(max(real(PEz))) max(max(real(PEz)))]);colorbar;
%     figure(); pcolor(zz_mesh,rr_mesh,(real(PHphi))); shading interp; hold on; contour(zz_mesh,rr_mesh,(real(PHphi))); hold off;
%     xlabel('Length [m]'); ylabel('Radius [m]'); title('Field H_{phi} [A/m]');caxis([-max(max(real(PHphi))) max(max(real(PHphi)))]);colorbar;
    

    % Write fields
    
    nome_fileEz=['MMMlong_Ez_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];
    nome_fileEr=['MMMlong_Er_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];
    nome_fileHphi=['MMMlong_Hphi_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];
    nome_fileMeshGridZ=['MMMlong_MeshGridZ_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];
    nome_fileMeshGridR=['MMMlong_MeshGridR_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_Material_',material.name,...
             '_f',num2str(f),'_P',num2str(P),'_S',num2str(S),'.txt'];   
    dlmwrite([SaveDir,nome_fileEz],PEz,'precision',10);
    dlmwrite([SaveDir,nome_fileEr],PEr,'precision',10);
    dlmwrite([SaveDir,nome_fileHphi],PHphi,'precision',10);
    dlmwrite([SaveDir,nome_fileMeshGridZ],[zz_mesh],'precision',10);
    dlmwrite([SaveDir,nome_fileMeshGridR],[rr_mesh],'precision',10);
    disp([nome_fileEz]);
    disp([nome_fileEr]);
    disp([nome_fileHphi]);
    disp([nome_fileMeshGridR]);
    disp([nome_fileMeshGridZ]);
    end
end

%% Save
%  Define files
nome_fileRe=['MMMlong_Re_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_Material_',material.name,...
         '_fmin',num2str(min(material.freq)),...
         '_fmax',num2str(max(material.freq)),'_P',num2str(P),'_S',num2str(S),'.txt'];
nome_fileIm=['MMMlong_Im_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t',num2str(t),'_Material_',material.name,...
         '_fmin',num2str(min(material.freq)),...
         '_fmax',num2str(max(material.freq)),'_P',num2str(P),'_S',num2str(S),'.txt'];
nome_filePP=['MMMlong_Postproc_L',num2str(L),'_Beta',num2str(beta),'_b',num2str(b),'_t','_Material_',material.name,...
         '_fmin',num2str(min(material.freq)),...
         '_fmax',num2str(max(material.freq)),'_P',num2str(P),'_S',num2str(S),'.txt'];
  
% Write impedance
dlmwrite([SaveDir,nome_fileRe],[material.freq',real(Z)'],'precision',10);
dlmwrite([SaveDir,nome_fileIm],[material.freq',imag(Z)'],'precision',10);
% Write Power Loss
if loss==1
    fileID=fopen([SaveDir,nome_filePP],'w');
    for ii=1:length(f_post)
        fprintf(fileID,'%s %f\n%s %6.2f\n%s %6.2f\n%s %6.2f\n%s %6.2f\n%s %6.2f\n%s %6.2f\n%s %6.2f\n%s %6.2f\n',...
            'Frequency_[GHz]:', f_post(ii)/1e9,...
            'Real(Poynting)_[W]:',real(W1(ii)),...
            'Imag(Poynting)_[W]:',imag(W1(ii)),...
            'Conductivity_Load_losses_[W]:',real(W2(ii)),...
            'Real(Electric_Load_Losses)_[W]:',real(W4(ii)),...
            'Imag(Electric_Load_Losses)_[W]:',imag(W4(ii)),...
            'Real(Magnetic_Load_Losses)_[W]:',real(W3(ii)),...
            'Imag(Magnetic_Load_Losses)_[W]:',imag(W3(ii)),...
            'Real(Electric_Cavity_Losses)_[W]:',real(W7(ii)),...
            'Imag(Electric_Cavity_Losses)_[W]:',imag(W7(ii)),...
            'Real(Magnetic_Cavity_Losses)_[W]:',real(W6(ii)),...
            'Imag(Magnetic_Cavity_Losses)_[W]:',imag(W6(ii)),...
            'Q_total: ',Qvec(ii));
    end
    fclose(fileID);
end
disp(['Impedance data saved in: ',SaveDir]);
disp([nome_fileRe]);
disp([nome_fileIm]);
if loss==1
disp(nome_filePP)  
end
end % function ends