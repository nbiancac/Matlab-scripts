function [fres,beta,q0,loss,unc_fres,unc_beta,unc_q0,unc_loss]=S11fit_attenuation(freq11,xx11,yy11,startbeta,attenuation)
% Fitting of the reflection parameter for a single resonance.
% Fitting the data as a function of the normalised frequency shift delta
% (delta = freq/fres-fres/freq) to find beta, q0 (as in the S11fit_delta 
% routine). The resonance frequency is
% then determined by fitting the whole S11 formula, using the previous
% value as startpoint. Matlab curve fitting toolbox is needed.
%
% free parameters of the fit: resonant frequency, beta, Q0, attenuation
% fres      resonant frequency
% beta      coupling factor
% q0        unloaded quality factor
% loss      constant attenuation of the feeding trasmission line
% unc_fres  uncertainty on the resonant frequency
% unc_beta  uncertainty on the coupling factor
% unc_q0    uncertainty on the unloaded quality factor
% unc_loss  uncertainty on constant attenuation of the feeding trasmission line
% Input parameters
% freq11    frequency (independent variable)
% xx11      real part of the S11 to be fitted
% yy11      imag part of the S11 to be fitted
% startbeta first guess of the copling factor <1 (>1) for under (over) coupling
% attenuation   constant attenuation of the feeding trasmission line
%               (fixed parameter in the fit)
% function [fres,beta,q0,loss,unc_fres,unc_beta,unc_q0,unc_loss]=S11fit_attenuation(freq11,xx11,yy11,startbeta,attenuation)
% ver. 1.1 Andrea Mostacci 07/10/08 

% to allow column and raw vector as input 
[raw,col]=size(freq11);
if col>raw,
    freq11=freq11';
end

[raw,col]=size(xx11);
if col>raw,
    xx11=xx11';
end

[raw,col]=size(yy11);
if col>raw,
    yy11=yy11';
end


mod_S11=sqrt(xx11.^2+yy11.^2);
[val,ind1]=min(mod_S11);
f0_S11=freq11(ind1);
delta_S11=freq11./f0_S11-f0_S11./freq11;

% Fitting as function of the normalised frequency shift
myS11 = fittype('factor*sqrt((((beta-1)/(beta+1))^2+(x*q0/(1+beta))^2)/(1+(x*q0/(1+beta))^2))',...
    'independent',{'x'},'coefficients',{'beta','q0','factor'});

% Fitting parameters for the Under/Over coupled cases
if startbeta <= 1,
    startbeta=(1-mod_S11(ind1))/(1+mod_S11(ind1));
    fitParam = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[0 0 0],'Upper',[1 Inf 1],'Startpoint',[startbeta 10000 attenuation]);
else
    startbeta=(1+mod_S11(ind1))/(1-mod_S11(ind1));
    fitParam = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[1 0 0],'Upper',[Inf Inf 1],'Startpoint',[startbeta 10000 attenuation]);
end

% Results of the first fit as a function of the normalised frequency
[fitS11,gofS11,outS11] = fit(delta_S11,mod_S11,myS11,fitParam);

coeffS11=coeffvalues(fitS11);
beta=coeffS11(1);
q0=coeffS11(2);
loss=coeffS11(3);

% Finding the resonance frequency using the previos results as start points
myS11freq = fittype('factor*sqrt((((beta-1)/(beta+1))^2+((x/f0-f0/x)*q0/(1+beta))^2)/(1+((x/f0-f0/x)*q0/(1+beta))^2))',...
    'independent',{'x'},...
             'coefficients',{'beta','q0','f0','factor'});
Fit4Freq = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[beta*0.9 q0*0.9 min(freq11) 0.5*loss],'Upper',[beta*1.1 q0*1.1 max(freq11) 1.5*loss],'Startpoint',[beta q0 f0_S11 loss]);


[fitFreq,gofFreq,outFreq] = fit(freq11,mod_S11,myS11freq,Fit4Freq);
coeffS11Freq=coeffvalues(fitFreq);
deltaParametersFreq=confint(fitFreq,0.95);

[aux,Nparameters]=size(deltaParametersFreq);
Ndegreees_of_freedom=length(freq11)-Nparameters;
t_value=tinv(0.95,Ndegreees_of_freedom);


beta=coeffS11Freq(1);
q0=coeffS11Freq(2);
fres=coeffS11Freq(3);
loss=coeffS11Freq(4);
unc_beta=(deltaParametersFreq(2,1)-deltaParametersFreq(1,1))./t_value;
unc_q0=(deltaParametersFreq(2,2)-deltaParametersFreq(1,2))./t_value;
unc_fres=(deltaParametersFreq(2,3)-deltaParametersFreq(1,3))./t_value;
unc_loss=(deltaParametersFreq(2,4)-deltaParametersFreq(1,4))./t_value;
