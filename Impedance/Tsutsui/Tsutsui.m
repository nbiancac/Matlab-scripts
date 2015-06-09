function [Zlong,Zxdriv,Zxdet,Zydriv,Zydet]=Tsutsui(a,b,d,L,freqN,material)
%
%Rectangular model for kicker using Tsutsui formalism
%
%function [Zlong,Zxdriv,Zxdet,Zydriv,Zydet]=Tsutsui(a,b,d,L,freqN)
%a=half width [m]
%b=half heigth (center-ferrite border) [m]
%d=half heigth (center-PEC border including ferrite) [m]
%L=Kicker's length [m]
%



u0=4*pi*1e-7;
material.rho=1./material.sigma;
e0=8.84*(1e-12);
c=299792458;
Z0=c*u0;
n=[0:1:20];
nx=[1:1:20];
Kxnx=(((nx))*pi)/(a);
Kxn=((2*n+1))*pi/(2*a);
sh=sinh(Kxn*b);
shx=sinh(Kxnx*b);
ch=cosh(Kxn*b);
chx=cosh(Kxnx*b);

Zlong=[];
Zxdriv=[];
Zydriv=[];
Zxdet=[];
Zydet=[];

if size(freqN,2)>size(freqN,1) 
    freqN=freqN';
end

for jj=1:1:length(freqN)
f=freqN(jj,1)
K=(2*pi*f)/c;
%  ur=1+(ui./(1+(k*k.*f.*f)))-1i*((ui*k.*f)./(1+(k*k.*f.*f)));
ur=material.mu_r(jj);
er=material.e_r(jj)-1i*(1./(2*pi*f*material.rho(jj)*e0)); %ferrite
%er=1-1i*1e5./(2*pi*f*e0); ur=1;   %graphite

Kyn=((((er.*ur)-1).*(K.^2))-(Kxn.^2)).^0.5;
Kynx=((((er.*ur)-1).*(K.^2))-(Kxnx.^2)).^0.5;
tn=tan(Kyn*(b-d));
tnx=tan(Kynx*(b-d));
ct=cot(Kyn*(b-d));
ctx=cot(Kynx*(b-d));
FX=(Kxn./K).*(1+(er.*ur)).*sh.*ch./((er.*ur)-1);
FXx=(Kxnx./K).*(1+(er.*ur)).*shx.*chx./((er.*ur)-1);
FY=(Kyn./K).*((ur.*(sh.^2).*tn)-(er.*(ch.^2).*ct))./((er.*ur)-1);
FYx=(Kynx./K).*((ur.*(shx.^2).*tnx)-(er.*(chx.^2).*ctx))./((er.*ur)-1);
FYv=(Kyn./K).*((ur.*(ch.^2).*tn)-(er.*(sh.^2).*ct))./((er.*ur)-1);
ImpLN=L*1i*(Z0/(2*a)).*sum((FX+FY-((K./Kxn).*sh.*ch)).^-1);
ImpHN=L*1i*(Z0/(2*a)).*sum(((Kxnx.^2)./K).*((FXx+FYx-((K./Kxnx).*shx.*chx)).^-1));
ImpVN=L*1i*(Z0/(2*a)).*sum(((Kxn.^2)./K).*((FX+FYv-((K./Kxn).*sh.*ch)).^-1));
ImpHdetN=-L*1i*(Z0/(2*a)).*sum(((Kxn.^2)./K).*(FX+FY-((K./Kxn).*sh.*ch)).^-1);
ImpVdetN=L*1i*(Z0/(2*a)).*sum(((Kxn.^2)./K).*(FX+FY-((K./Kxn).*sh.*ch)).^-1);
Zlong=[Zlong ImpLN];
Zxdriv=[Zxdriv ImpHN];
Zydriv=[Zydriv ImpVN];
Zxdet=[Zxdet ImpHdetN];
Zydet=[Zydet ImpVdetN];
end
