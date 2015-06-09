%
% Plot bessel functions
%
close all;
t=0:2/1000:10;

figure(1);
plot(t,besselk(1,t),'b','linewidth',2);
hold on;
plot(t,1./t,'r','linewidth',2);
%plot(t,besselk(2,t),'g','linewidth',1);
%plot(t,besselk(3,t),'m','linewidth',1);
grid on;
%legend('K0','K1','K2','K3');
%ylim([0 100])
%xlim([0 0.01])
figure
err=abs(besselk(1,t)-1./t)*100./besselk(1,t);
plot(t,err,'b','linewidth',2);
ylim([0 50]);
grid on;
title('K_1(z) error z\rightarrow0','fontsize',14)
ylabel('% error','fontsize',14);
xlabel('z','fontsize',14);
%% sen cos sinh
t=0:pi/100:pi/2-pi/2/100;
figure(3)
plot(t,(cos(t)-(1))*100./cos(t),'b','linewidth',2);

%%
WangDir='/afs/cern.ch/user/n/nbiancac/scratch0/HEADTAIL_RELEASE/Wang';
cd(WangDir)

Z0=120*pi;
q=1.602176e-19;
c=299792458;
u0=Z0/c;
e0=8.84*(1e-12);
gamma=450; %SPS inj
beta=sqrt(1-1/gamma^2)
v=beta*c;
a=0.05; %[m]
b=0.0305; %[m]
d=b+0.06; %[m]
L=2.738; %[m]
Np=100; % Numero di punti

x=(0.0001:(a-0.0001)/Np:a);
y=(0.0001:(b-0.0001)/Np:b);
q=1.6e-19;
f=10e6; %[Hz]
k0=2*pi.*f/c
kr=2*pi.*f/(v*gamma)

ExS=zeros(Np+1,Np+1);
EyS=zeros(Np+1,Np+1);
EzS=zeros(Np+1,Np+1);
C=0;

for i=1:1:length(x);
    for j=1:1:length(y);
        %r=sqrt((x(i)-2*n*a)^2+y(j)^2);
        %ExS(i,j)=Z0*q/2/pi/beta*besselk(1,kr*r)*kr*x(i)/r;
        %EyS(i,j)=EyS(i,j)+(-1)^n*Z0*q/2/pi/beta*kr*besselk(1,kr*r)*y(j)/r;
        EyS(i,j)=Z0*q/2/pi/beta*2*pi/2/a*(cos(pi*x(i)/2/a)*sinh(pi*y(j)/2/a))/(cosh(pi*y(j)/a)-cos(pi*x(i)/a));
        %EzS(i,j)=Z0*q/2/pi/beta/gamma*besselk(0,kr*r)*kr;
    end
end


figure(3)
surface(x,y,EyS); shading 'interp'

n=[0:1:20];
Zlong=[];
freq=[];
rho=1e6;
ui=460;
kf=1/(20e6);

for f=1e7:1e7:100e9
k0=2*pi*f/c;
k=2*pi*f/beta/c;
kr=2*pi*f/(beta*c*gamma);

N=Z0*q/2/beta/a;
er=12-1i*(1./(2*pi*f*rho*e0));
ur=1+(ui./(1+(kf*kf.*f.*f)))-1i*((ui*kf.*f)./(1+(kf*kf.*f.*f)));

g1n=(2*n+1)*pi/2/a;
k1n=(2*n+1)*pi/2/a;
k2n=sqrt(-kr^2-k1n.^2);
g2n=sqrt(k0^2*ur*er-k^2-g1n.^2);
x1=k2n*b;
x2=g2n*(b-d);

An = N*exp(1i.*x1).*( k2n.*cos(x1).*(1i*k2n.*gamma^2*(k^2-k0^2*er*ur)-g2n.*k^2*er.*cot(x2)) + sin(x1).*...
    (k^2.*(kr^2*er*ur+k1n.^2*(1+gamma^2*(-1+er*ur))) +1i.*g2n.*k2n*kr^2*gamma^2*ur.*tan(x2)) )./...
    (k.*k2n*gamma^2.*(-g2n.*k2n*er.*cos(x1).^2.*cot(x2)+(k^2+(-k0^2+kr^2)*er*ur+k1n.^2*(1+er*ur)).*...
    cos(x1).*sin(x1)-g2n.*k2n*ur.*sin(x1).^2.*tan(x2)));

Zlong=[Zlong, -L/q*sum(An)];
freq=[freq,f];
end

figure(4)
plot(freq./1e9,real(Zlong),'b','LineWidth',2);
hold on
plot(freq./1e9,imag(Zlong),'r','LineWidth',2)
plot(ftsu,real(Zlongtsu),'--r','LineWidth',2)
plot(ftsu,imag(Zlongtsu),'--b','LineWidth',2)
legend('Re(Zs) Wang',   'Im(Zs) Wang',...
       'Re(Zs) Tsutsui','Im(Zs) Tsutsui');
xlim([0 4]);
title(['Tsutsui-Wang Longitudinal Impedance Zs for \beta=1'],'fontsize',14)
xlabel('Frequency in GHz','fontsize',14)
ylabel('Impedance in \Omega /m','fontsize',14)
grid on

saveas(gcf,['Tsutsui-Wang.fig'],'fig');
saveas(gcf,['Tsutsui-Wang.pdf'],'pdf');

%%
% Serie approssimata
a=0.05;
x=0.04;
y=0.01;
C=0
S1=y/(x^2+y^2);
n=1;
err=100;
while err>1e-10

S1=S1+(-1)^n*(y)/((x-2*n*a)^2+y^2)+(-1)^(-n)*(y)/((x-2*(-n)*a)^2+y^2);
C=[C;S1];
err=abs(C(end)-C(end-1));
n=n+1;
end
S1
n

% Serie Bessel
C2=0
r=sqrt((x)^2+y^2);
S3=besselk(1,r)*y/r
n=1;
err=100;
while err>1e-10

S3=S3+(-1)^n*(y)/sqrt((x-2*n*a)^2+y^2)*besselk(1,sqrt((x-2*n*a)^2+y^2))+(-1)^-n*(y)/sqrt((x-2*n*a)^2+y^2)*besselk(1,sqrt((x+2*n*a)^2+y^2));
C2=[C2;S3];
err=abs(C2(end)-C2(end-1));
n=n+1;
end
S3
n

% Serie Somma
S2=2*pi/2/a*(cos(pi*x/2/a)*sinh(pi*y/2/a))/(cosh(pi*y/a)-cos(pi*x/a))

%%
a=0.05;
b=0.002;
x=(-a:2*a/1000:+a);
I=0;
%trapz( x , cos(pi*x./(2*a))*sinh(pi*b/(2*a))  /  (   cosh(pi*b/a) - cos(pi*x/a)    ) * cos((2*n+1)*pi.*x/2/a)  )
 b=0.000000012
yy=pi*b/2/a./(0.5*(pi*b/a)^2+2*sin(pi.*x/2/a).^2);
plot(x,yy)
I=[I; trapz(x,yy)];
