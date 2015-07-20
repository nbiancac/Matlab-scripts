function LHC=LHC_param(E0,E,scenario)
% parameters for the LHC
% 04-06-2015: the bunch length has been changed to 9cm from 7.5cm as delivered from the injectors.        


LHC={};
LHC.chromay=0;
LHC.chromax=0;
c=299792458;
if E==7e12
    if strcmp(scenario,'Nominal LHC')
        LHC.Estr='7TeV';
        LHC.Nb=1.15e11;
        LHC.M=2808;
        LHC.sigmaz=9e-2; % rms
        LHC.V=16e6;
        LHC.scenario={scenario};
    end
elseif E==6.5e12
    if strcmp(scenario,'Nominal LHC')
        LHC.Estr='6500GeV';
        LHC.Nb=1.15e11;
        LHC.M=2808;
        LHC.sigmaz=9e-2; % rms
        LHC.V=16e6;
        LHC.scenario={scenario};
    elseif strcmp(scenario, 'RF-MD-12MV')
        LHC.Estr='6500GeV';
        LHC.Nb=1.6e11;
        LHC.M=1;
        LHC.sigmaz=6e-2; % rms
        LHC.V=12e6;
        LHC.scenario={scenario};
    elseif strcmp(scenario, 'RF-MD-6MV')
        LHC.Estr='6500GeV';
        LHC.Nb=1.6e11;
        LHC.M=1;
        LHC.sigmaz=7.12e-2; % rms
        LHC.V=6e6;
        LHC.scenario={scenario};
    end
elseif E==4e12
    if strcmp(scenario,'Nominal LHC')
        LHC.Estr='4000GeV';
        LHC.Nb=1.15e11;
        LHC.M=2808;
        LHC.sigmaz=1.25*9e-2; % rms
        LHC.V=12e6;
        LHC.scenario={scenario};
    end
elseif E==450e9
    if strcmp(scenario,'Nominal LHC')
        LHC.Estr='450GeV';
        LHC.Nb=1.15e11*1.05;
        LHC.M=2808;
        LHC.sigmaz=1.3*9e-2; % rms
        LHC.V=6e6;
        LHC.scenario={scenario};
    end    
end
    

e=1.602176487e-19; 
c=299792458;

LHC.machine='LHC';
LHC.h=35640; 
LHC.gamma=E*e/E0;
LHC.beta=sqrt(1.-1./(LHC.gamma^2));
LHC.mom=LHC.gamma*LHC.beta*E0/c/(e/c); % momentum eV/c
LHC.taub=LHC.sigmaz/LHC.beta/c; % rms bunch length
LHC.circ=26658.883; % total circumference in m
LHC.R=LHC.circ/(2.*pi); % machine radius
LHC.Qx=64.31;LHC.Qxfrac=LHC.Qx-floor(LHC.Qx);
LHC.Qy=59.32;LHC.Qyfrac=LHC.Qy-floor(LHC.Qy);
LHC.alphap=3.225e-4; % momentum compaction factor
LHC.eta=LHC.alphap-1./(LHC.gamma^2); % slip factor
LHC.gammatr=sqrt(abs(1/LHC.alphap));
LHC.phis=0;
LHC.Qs=Qs_from_RF_param(LHC.V,LHC.h,LHC.gamma,LHC.eta,LHC.phis,'proton');
LHC.f0=c*LHC.beta/LHC.circ; % rev. frequency
LHC.omega0=2.*pi*LHC.f0;
LHC.omegas=LHC.Qs*LHC.omega0;
LHC.dphase=0.; 
LHC.betavx=LHC.R/LHC.Qx;
LHC.betavy=LHC.R/LHC.Qy;
    
end

