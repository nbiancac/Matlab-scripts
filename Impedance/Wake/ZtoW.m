function ZtoW(name,flag_show)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   description: calculates wake from impedance with ifft.
%   input: name of the .imp file
%   output: name.wake file and pictures in the Wakes directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% directories
mainDir='/afs/cern.ch/project/impedance/PS_impedance_database/';
ImpDir=[mainDir,'Impedances/'];
WakeDir=[mainDir,'Wakes/'];
CodesDir=[mainDir,'Codes/'];
ResultsDir=[mainDir,'Results/'];
OpticsDir=[mainDir,'Optics/'];
BeamDir=[mainDir,'Beams/']; 
path(path,CodesDir);

% constants
c=299792458; 

% impedance
disp('Transverse wake calculation...')
if exist([ImpDir,'Transverse/',name,'.imp'],'file')
    L=dlmread([ImpDir,'Transverse/',name,'.imp'],'',1,0);
    f=L(:,1)';
    Zx=L(:,2)'+1j*L(:,3)';
    Zy=L(:,4)'+1j*L(:,5)';
    Zxdet=L(:,6)'+1j*L(:,7)';
    Zydet=L(:,8)'+1j*L(:,9)';

    % double the spectrum
    if f(1)>0
        fforFFT=[-f(end:-1:1) 0 f];
        ZydetforFFT=[-conj(Zydet(end:-1:1))  1i*imag(Zydet(1)) Zydet];
        ZyforFFT=[-conj(Zy(end:-1:1))  1i*imag(Zy(1)) Zy];
        ZxdetforFFT=[-conj(Zxdet(end:-1:1))  1i*imag(Zxdet(1)) Zxdet];
        ZxforFFT=[-conj(Zx(end:-1:1))  1i*imag(Zx(1)) Zx];
    elseif f(1)==0
        fforFFT=[-f(end:-1:2) f];
        ZydetforFFT=[-conj(Zydet(end:-1:2)) Zydet];
        ZyforFFT=[-conj(Zy(end:-1:2)) Zy];
        ZxdetforFFT=[-conj(Zxdet(end:-1:2))  Zxdet];
        ZxforFFT=[-conj(Zx(end:-1:2)) Zx];
    end
    
    % interpolation
    if std(diff(f))==0 
        disp('Frequencies  are equispaced.')
        fint=f;
    else
        fstep=1e6;
        fint=[min(f):fstep:max(f)];
        disp('Frequencies not equispaced: interpolating...')
        disp(['Step frequency is ',num2str(min(fstep/1e6)),' MHz']);
        disp(['Min frequency is ',num2str(min(f/1e6)),' MHz']);
        disp(['Max frequency is ',num2str(max(f/1e9)),' GHz']);
    end
    
    fint=[-fint(end:-1:1) 0 fint];
    ZxforFFTint=interp1(fforFFT,ZxforFFT,fint);
    ZyforFFTint=interp1(fforFFT,ZyforFFT,fint);
    ZxdetforFFTint=interp1(fforFFT,ZxdetforFFT,fint);
    ZydetforFFTint=interp1(fforFFT,ZydetforFFT,fint);
    
    % Ifft
    N=size(fint,2);
    df=(fint(2)-fint(1));
    dt=1/(N*df);% sec
    time=(0:1:N-1)*dt;
    time=time-max(time)/2;
    Wydetk=[time;-df*N*imag(fftshift(ifft(ifftshift(ZydetforFFTint))))]';    % (s,V/Cm)
    Wyk=   [time;-df*N*imag(fftshift(ifft(ifftshift(ZyforFFTint))))]';          % (s,V/Cm)
    Wxdetk=[time;-df*N*imag(fftshift(ifft(ifftshift(ZxdetforFFTint))))]';    % (s,V/Cm)
    Wxk=   [time;-df*N*imag(fftshift(ifft(ifftshift(ZxforFFTint))))]';          % (s,V/Cm)
    
    % prepare for Headtail (N.B. wake sign change)
    WyHDTLk=   [Wyk(:,1)    -Wyk(:,2)*1e-12/1000]; % (s,V/pCmm)
    WxHDTLk=   [Wxk(:,1)    -Wxk(:,2)*1e-12/1000]; % (s,V/pCmm)
    WydetHDTLk=[Wydetk(:,1) -Wydetk(:,2)*1e-12/1000]; % (s,V/pCmm)
    WxdetHDTLk=[Wxdetk(:,1) -Wxdetk(:,2)*1e-12/1000]; % (s,V/pCmm)

    figure();
    set(gcf,'visible',flag_show)
    plot(time,(WxHDTLk(:,2)),'-k'); hold on;
    plot(time,(WxdetHDTLk(:,2)),'--k'); 
    plot(time,(WyHDTLk(:,2)),'-r');
    plot(time,(WydetHDTLk(:,2)),'--r'); hold off;
    legend('W_x^{dip}','W_x^{det}','W_y^{dip}','W_y^{det}');
    xlim([-5 20].*1e-10)
    xlabel('time [s]')
    ylabel('W_t [V/pCmm]')
    s=load([CodesDir,'Plots/','style.mat']);;
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [WakeDir,'Transverse/',name,'.fig'],'fig');
    hgexport(gcf, [WakeDir,'Transverse/',name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [WakeDir,'Transverse/',name,'.png'],s,'Format','png');
    
    % Writing 4 wake for HDTLattice
    wakes4HDTLk=[WxHDTLk(:,1),WxHDTLk(:,2),WyHDTLk(:,2),WxdetHDTLk(:,2),WydetHDTLk(:,2)];
    dlmwrite([WakeDir,'Transverse/',name,'.wake'],wakes4HDTLk,'delimiter','\t');
    disp('Wakes written in Transverse:')
    disp([name,'.wake'])
else
    warning([name,'.imp',' file do not exist in ',ImpDir]);
    warning('Transverse wake calculation aborted.')
end % end transverse

disp('Longitudinal wake calculation...')
if exist([ImpDir,'Longitudinal/',name,'.imp'],'file')
    L=dlmread([ImpDir,'Longitudinal/',name,'.imp'],'',1,0);
    f=L(:,1)';
    Zlong=L(:,2)'+1j*L(:,3)'; 
    
    % double the spectrum
    fforFFT=[-f(end:-1:1) 0 f];
    ZlongforFFT=[-conj(Zlong(end:-1:1))  1i*imag(Zlong(1)) Zlong];
    
    if std(diff(f))==0 
        disp('Frequencies  are equispaced.')
        fint=f;
    else
        fstep=1e7;
        fint=[min(f):fstep:max(f)];
        disp('Frequencies not equispaced: interpolating...')
        disp(['Step frequency is ',num2str(min(fstep/1e6)),' MHz']);
        disp(['Min frequency is ',num2str(min(f/1e6)),' MHz']);
        disp(['Max frequency is ',num2str(max(f/1e9)),' GHz']);
    end
    fint=[-fint(end:-1:1) 0 fint];
    ZlongforFFTint=interp1(fforFFT,ZlongforFFT,fint);

    % Ifft
    N=size(fint,2);
    df=(fint(2)-fint(1));
    dt=1/(N*df);% sec
    time=(0:1:N-1)*dt;
    time=time-max(time)/2;
    Wlongk=[time;-df*N*imag(fftshift(ifft(ifftshift(ZlongforFFTint))))]'; % (s,V/C)

    % prepare for Headtail (N.B. Wake with sign change)
    WlongHDTLk=[Wlongk(:,1) -Wlongk(:,2)*1e-12]; % (s, V/pC)

    figure(2);
    set(gcf,'visible',flag_show)
    plot(time,(WlongHDTLk(:,2)),'-k'); 
    legend('W_{long}');
    xlim([-5 20].*1e-10)
    xlabel('time [s]');
    ylabel('W_l [V/pC]')
    s=load([CodesDir,'Plots/','style.mat']);
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [WakeDir,'Longitudinal/',name,'.fig'],'fig');
    hgexport(gcf, [WakeDir,'Longitudinal/',name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [WakeDir,'Longitudinal/',name,'.png'],s,'Format','png');
        
    wakes4HDTLk=[WlongHDTLk(:,1),WlongHDTLk(:,2)];
    dlmwrite([WakeDir,'Longitudinal/',name,'.wake'],wakes4HDTLk,'delimiter','\t');
    disp('Wakes written in Longitudinal:')
    disp([name,'.wake'])
else
    warning([name,'.imp',' file do not exist in ',ImpDir]);
    warning('Longitudinal wake calculation aborted.')
end % end longitudinal

end % end function
