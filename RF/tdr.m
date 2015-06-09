function [pulse_out,step_out,Z_step]=tdr(Z0,sigma_t,MeasDir,ResultDir,name,saveplot)
%   function [pulse_out,step_out,Z_step]=tdr(Z0,sigma_t,MeasDir,ResultDir,name,saveplot)
%   
%
    
    letto=importdata([MeasDir,'TS0.S2P']);
    freq=letto.data(:,1);
    S11=(letto.data(:,2).*exp(1i.*letto.data(:,3)*pi/180));
    S21=(letto.data(:,4).*exp(1i.*letto.data(:,5)*pi/180));
    S12=(letto.data(:,6).*exp(1i.*letto.data(:,7)*pi/180));
    S22=(letto.data(:,8).*exp(1i.*letto.data(:,9)*pi/180));

%     figure(1)
%     plot(freq,20*log10(abs(S11)));
%     xlabel('f');
%     ylabel('S11 dB')

    N=numel(freq);
    freq_int=linspace(0,max(freq),N);
    S11_int=interp1(freq,S11,freq_int,'linear','extrap');
    S22_int=interp1(freq,S22,freq_int,'linear','extrap');
    S12_int=interp1(freq,S12,freq_int,'linear','extrap');
    S21_int=interp1(freq,S21,freq_int,'linear','extrap');

    Nfull=2*N-1;
    freq_full=[-freq_int(end:-1:2) freq_int(1:end)];
    S11_full=[conj(S11_int(end:-1:2)) real(S11_int(1)) S11_int(2:end)];
    S22_full=[conj(S22_int(end:-1:2)) real(S22_int(1)) S22_int(2:end)];
    S12_full=[conj(S12_int(end:-1:2)) real(S12_int(1)) S12_int(2:end)];
    S21_full=[conj(S21_int(end:-1:2)) real(S21_int(1)) S21_int(2:end)];
    
    df=diff(freq_full(1:2));
    dt=1/(df*Nfull);
    time_pulse=dt.*[(-Nfull+1)/2:1:(Nfull-1)/2];
    clight=299792456.2;
    space_pulse=clight.*time_pulse;
    
    wind=1./sqrt(2*pi)/sigma_t*exp(-1/2/sigma_t^2*(time_pulse).^2);
    wind_freq=exp(-0.5*(2*pi*freq_full.*sigma_t).^2);
    
    S11t_pulse=real(ifftshift(ifft(fftshift(S11_full.*wind_freq))));
    S22t_pulse=real(ifftshift(ifft(fftshift(S22_full.*wind_freq))));
    S12t_pulse=real(ifftshift(ifft(fftshift(S12_full.*wind_freq))));
    S21t_pulse=real(ifftshift(ifft(fftshift(S21_full.*wind_freq))));
    
    step=ones(Nfull,1);
    space_step=[space_pulse,2*space_pulse(end)+space_pulse(2:end)];
    S11t_step=conv(step,S11t_pulse);
    S22t_step=conv(step,S22t_pulse);
    S12t_step=conv(step,S12t_pulse);
    S21t_step=conv(step,S21t_pulse);

    if saveplot==1
        figure(2);
        plot(space_step,S11t_step); hold on;
        plot(space_step,S22t_step,'r'); hold on;
        legend('S_{1,1}','S_{2,2}')
        xlabel('s [m]')
        ylabel('S(s)')
        xlim([-1 8])
        name=[name,'_step'];
        s=hgexport('readstyle','PRSTAB');
        hgexport(gcf,'',s,'applystyle',true);
        saveas(gcf, [ResultDir,name,'.fig'],'fig');
        hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
        hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
    end
    
    pulse_out=[space_pulse',S11t_pulse',S21t_pulse',S12t_pulse',S22t_pulse'];
    step_out=[space_step',S11t_step',S21t_step',S12t_step',S22t_step'];
    
    Z_step=zeros(size(step_out));
    Z_step(:,1)=step_out(:,1);
    Z_step(:,2)=Z0*(1+step_out(:,2))./(1-step_out(:,2));
    Z_step(:,3)=Z0*(1+step_out(:,3))./(1-step_out(:,3));
    Z_step(:,4)=Z0*(1+step_out(:,4))./(1-step_out(:,4));
    Z_step(:,5)=Z0*(1+step_out(:,5))./(1-step_out(:,5));
    
end