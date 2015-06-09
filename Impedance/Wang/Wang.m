%% Nang Wang formula for Kickers with beta<1
%
%   Wang(0.05,0.0305,0.0305+0.06,2.73,27.7, 1e6:100e6:1e9)
%
%   Rectangular model for kicker using Tsutsui-Wang formalism
%
%   a=half width [m]
%   b=half heigth (center-ferrite border) [m]
%   d=half heigth (center-PEC border including ferrite) [m]
%   L=Kicker's length [m]
%   gamma=relativistic gamma
%   freq: frequencies in which calculate impedances
%


function [Zlong,Zxdip,ZquadX,Zydip,ZquadY]=Wang(a,b,d,L,gamma, freq, material)

Zlong=[];
Zxdip=[];
ZquadX=[];
Zydip=[];
ZquadY=[];

%% Physical constants and material data

Z0=120*pi;
q=1.602176e-19;
c=299792458;
u0=Z0/c;
e0=8.84*(1e-12);
beta=sqrt(1-1/gamma^2)
v=beta*c;

if strcmp(material,'ferrite')
    rho=1e6;
    ui=460;
    kf=1/(20e6);
    media_ur='ur=1+(ui./(1+(kf*kf.*freq(jj).*freq(jj))))-1i*((ui*kf.*freq(jj))./(1+(kf*kf.*freq(jj).*freq(jj))));';
    media_er='er=12-1i*(1./(om*rho*e0));';
elseif strcmp(material,'graphite')
    media_ur='ur=1;';
    media_er='er=1-1i*1e5./(om*e0);';
end
    
%% Longitudinal + Quadrupolar impedance
Nmax=170;

N=Z0*q/2/beta/a;
An=zeros(Nmax);

for jj=1:1:length(freq)
    
n=0:Nmax;
g1n=(2*n+1)*pi/2/a;
k1n=(2*n+1)*pi/2/a;

om=2*pi*freq(jj);
k0=om/c;
k=om/beta/c;
kr=om/(beta*c*gamma);

eval(media_ur);
eval(media_er);

k2n=sqrt(-kr^2-k1n.^2);
g2n=sqrt(k0^2*ur*er-k^2-g1n.^2);
x1=k2n*b;
x2=g2n*(b-d);
Eyn=N*exp(1i.*x1);

An = Eyn.*( k2n.*cos(x1).*(1i*k2n.*gamma^2*(k^2-k0^2*er*ur)-g2n.*k^2*er.*cot(x2)) + sin(x1).*...
   (k^2.*(kr^2*er*ur+k1n.^2*(1+gamma^2*(-1+er*ur))) +1i.*g2n.*k2n*kr^2*gamma^2*ur.*tan(x2)) )./...
   (k.*k2n*gamma^2.*(-g2n.*k2n*er.*cos(x1).^2.*cot(x2)+(k^2+(-k0^2+kr^2)*er*ur+k1n.^2*(1+er*ur)).*...
   cos(x1).*sin(x1)-g2n.*k2n*ur.*sin(x1).^2.*tan(x2)));

Zlong=[Zlong;-L/q*sum(An)];

ZquadX=[ZquadX; L/k/q*sum(An.*k1n.^2)];
ZquadY=[ZquadY; L/k/q*sum(An.*k2n.^2)];
%u=0:1e-4:10;
%alpha00=2*trapz((1-exp(-2*k/gamma*a.*cosh(u)))./(2*sinh(2*k/gamma*a.*cosh(u))),u);
%Zsc=[Zsc,1i*om*Z0*L/(2*pi*beta^2*gamma^2*c)*alpha00];
end


%% Dipolar Impedance
% Y
Nmax=170;
YAn=zeros(Nmax);


for jj=1:1:length(freq)
    
n=0:Nmax;

k1n=(2*n+1)*pi/2/a;

om=2*pi*freq(jj);
k0=om/c;
k=om/beta/c;
kr=om/(beta*c*gamma);

eval(media_ur);
eval(media_er);

k2n=sqrt(-kr^2-k1n.^2);
k3n=sqrt(k0^2*ur*er-k^2-k1n.^2);
x1=k2n*b;
x2=k3n*(b-d);


M1n=1i*k2n.^2*k0.^2*(1-er*ur)+k2n.*kr^2.*(1i*k2n-k3n.*er.*cot(x2));
M2n=-k1n.^2*k0^2*(1-er*ur)+k2n.*kr^2*ur.*(-k2n.*er+1i*k3n.*tan(x2));
M3n=k2n.*k3n.*(er*sin(x1).^2.*cot(x2)+ur*cos(x1).^2.*tan(x2));
M4n=(k0^2*(-1+er*ur)+k2n.^2.*(1+er*ur)).*cos(x1).*sin(x1);

YAn=(1i*Z0*k2n*L/2/a/beta/k)  .*  k2n.*exp(1i.*k2n*b).*(-M1n.*sin(x1)+M2n.*cos(x1))./(k*k2n.*(M3n-M4n));

Zydip=[Zydip; sum(YAn)];

end

% X
Nmax=170;
XAn=zeros(Nmax);


for jj=1:1:length(freq)
    
n=0:Nmax;

k1n=(n)*pi/a;

om=2*pi*freq(jj);
k0=om/c;
k=om/beta/c;
kr=om/(beta*c*gamma);

eval(media_ur);
eval(media_er);

k2n=sqrt(-kr^2-k1n.^2);
k3n=sqrt(k0^2*ur*er-k^2-k1n.^2);
x1=k2n*b;
x2=k3n*(b-d);


M1n=1i*k2n.^2*k0.^2*(1-er*ur)+k2n.*kr^2.*(1i*k2n-k3n.*er.*cot(x2));
M2n=-k1n.^2*k0^2*(1-er*ur)+k2n.*kr^2*ur.*(-k2n.*er+1i*k3n.*tan(x2));
M3n=-k2n.*k3n.*(er*cos(x1).^2.*cot(x2)+ur*sin(x1).^2.*tan(x2));
M4n=(k0^2*(-1+er*ur)+k2n.^2.*(1+er*ur)).*cos(x1).*sin(x1);


XAn=(1i*Z0*L/2/a/beta/k)  *k1n.*exp(1i.*k2n*b).*(M1n.*cos(x1)+M2n.*sin(x1))./(k*(M3n-M4n)); 
%XAn=-L*(-1i*Z0*k2n)./(2*a*beta)  .*k1n.*  (-1i*k1n.^2*(1-er*ur))./(k2n.*(M3n-M4n));
Zxdip=[Zxdip; sum(XAn)];
end

