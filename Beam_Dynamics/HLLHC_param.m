function HLLHC=HLLHC_param(E0,E,scenario)

HLLHC={};
HLLHC.chromay=0;
HLLHC.chromax=0;

if E==7e12
    if strcmp(scenario,'HL-LHC 25ns')
        HLLHC.Estr='7TeV';
        HLLHC.Nb=2.2e11;
        HLLHC.M=2748;
        HLLHC.sigmaz=(1.0e-9/4)*constants('clight'); % rms
        HLLHC.V=16e6;
        HLLHC.scenario={scenario};
        HLLHC.en=2.5e-6;
        HLLHC.betavy=71.5255;
        HLLHC.betavx=65.9756;
    elseif strcmp(scenario,'HL-LHC BCMS')
        HLLHC.Estr='7TeV';
        HLLHC.Nb=2.2e11;
        HLLHC.M=2604;
        HLLHC.sigmaz=(1.0e-9/4)*constants('clight'); % rms
        HLLHC.V=16e6;
        HLLHC.en=2.5e-6;
        HLLHC.scenario={scenario};
    elseif strcmp(scenario,'HL-LHC 8b+4e')
        HLLHC.Estr='7TeV';
        HLLHC.Nb=2.2e11;
        HLLHC.M=1968;
        HLLHC.sigmaz=(1.0e-9/4)*constants('clight'); % rms        
        HLLHC.V=16e6;
        HLLHC.en=2.2e-6;
        HLLHC.scenario={scenario};
    end
elseif E==450e9
     if strcmp(scenario,'HL-LHC 25ns')
        HLLHC.Estr='450GeV';
        HLLHC.Nb=2.3e11;
        HLLHC.M=2748;
        HLLHC.sigmaz=(1.0e-9/4)*constants('clight'); % rms
        HLLHC.V=6e6;
        HLLHC.en=2; %micron normalized emittance
        HLLHC.scenario={scenario};
    elseif strcmp(scenario,'HL-LHC BCMS')
        HLLHC.Estr='450GeV';
        HLLHC.Nb=2.3e11;
        HLLHC.M=2604;
        HLLHC.sigmaz=(1.0e-9/4)*constants('clight'); % rms
        HLLHC.V=6e6;
        HLLHC.en=1.4; %micron normalized emittance
        HLLHC.scenario={scenario};
     elseif strcmp(scenario,'HL-LHC 8b+4e')
        HLLHC.Estr='450GeV';
        HLLHC.Nb=2.4e11;
        HLLHC.M=1968;
        HLLHC.sigmaz=(1.0e-9/4)*constants('clight'); % rms        
        HLLHC.V=6e6;
        HLLHC.en=1.7; %micron normalized emittance     
        HLLHC.scenario={scenario};
    end    
end
    

e=1.602176487e-19; 
c=constants('clight');

HLLHC.machine='HLLHC';
HLLHC.h=35640; 
HLLHC.gamma=E*e/E0;
HLLHC.beta=sqrt(1.-1./(HLLHC.gamma^2));
HLLHC.taub=HLLHC.sigmaz/HLLHC.beta/c; % rms bunch length
HLLHC.circ=26658.883; % total circumference in m
HLLHC.R=HLLHC.circ/(2.*pi); % machine radius
HLLHC.Qx=64.31;HLLHC.Qxfrac=HLLHC.Qx-floor(HLLHC.Qx);
HLLHC.Qy=59.32;HLLHC.Qyfrac=HLLHC.Qy-floor(HLLHC.Qy);
HLLHC.alphap=3.225e-4; % momentum compaction factor
HLLHC.eta=HLLHC.alphap-1./(HLLHC.gamma^2); % slip factor
HLLHC.gammatr=sqrt(abs(1/HLLHC.alphap));
HLLHC.phis=0;
HLLHC.Qs=Qs_from_RF_param(HLLHC.V,HLLHC.h,HLLHC.gamma,HLLHC.eta,HLLHC.phis,'proton');
HLLHC.f0=c*HLLHC.beta/HLLHC.circ; % rev. frequency
HLLHC.omega0=2.*pi*HLLHC.f0;
HLLHC.omegas=HLLHC.Qs*HLLHC.omega0;
HLLHC.dphase=0.; 
    
end

