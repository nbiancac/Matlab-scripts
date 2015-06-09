% Example2a 
%
% 2-Stub matching example 
% Using smdrawc.m, trl.m, and junct.m
%
% Load impedance Zload is matched using 2 stubs, 
% lengths LenS1 and LenS2 respectively. 
% Calculated over frequency range 1000 to 3000 MHz
%          
%
%         <-LenTx1->  |   <-LenTx2->   | <-LenTx3->
%     
% Zload =============Z1=Z1a===========Z2=Z2a======== Zin
%                Zstub1           Zstub2
%                     |                |
%                     |LenS1           |LenS2
%                     |                |
%                    S/C              S/C
%
%
% For this example : LenTx1=0 (mm)
%                    LenTx2=lambda/8 ( 18.75mm at 2000 MHz)
%                    LenTx3=0 (mm)
%                    (All tx-lines are 50 ohms, Er=1.0, Loss/m=0)
%                     
%
% Open this file and look through the comments.
%
% Ref. David M. Pozar 'Microwave Engineering 2nd Edition'   Worked Example 5.4

% N.Tucker www.activefrance.com 2008

clc;
close all;
fprintf('\n\n\n***** Example 2a *******\n\n\n');
help example2a;
% Set up some variables to use in the example
Zo=50;                           % 50 ohm impedance for smith chart
Freq=1000:10:3000;               % Frequency vector 1000 to 3000 MHz in 10MHz steps
Zload=term(60-80*j,Freq);        % Load impedance

MidFreq=(max(Freq)+min(Freq))/2; % Midpoint frequency (MHz)
MidLambda=3e8/(MidFreq*1e6)*1e3; % Midpoint wavelength (mm)

LenS1=0.232*MidLambda;           % Length of first stub (mm)
LenS2=0.100*MidLambda;           % Length of second stub (mm)

LenTx1=0;                        % Length of transmission line section 1
LenTx2=MidLambda/8;              % Length of transmission line section 2, between the stubs
LenTx3=0;                        % Length of transmission line section3

Er=1.0;                          % Dielectric const
LdB=0;                           % Loss in dB/m
ZoT=50;                          % Characteristic impedance (ohms)
Zsc=0;                           % Short circuit impedance (0 ohms)


% Model the various sections of tx-line and junctions

Zstub1=trl(ZoT,Zsc,LenS1,Freq,Er,LdB);   % Impedance at end of S/C stub-1
Zstub2=trl(ZoT,Zsc,LenS2,Freq,Er,LdB);   % Impedance at end of S/C stub-2

Z1=trl(ZoT,Zload,LenTx1,Freq,Er,LdB);    % Impedance on main tx-line at first juntion
Z1a=junct(Z1,Zstub1,Freq);               % Combined impedance of main-line and stub-1

Z2=trl(ZoT,Z1a,LenTx2,Freq,Er,LdB);      % Impedance on main tx-line at second juntion
Z2a=junct(Z2,Zstub2,Freq);               % Combined impedance of main-line and stub-2

Zin=trl(ZoT,Z2a,LenTx3,Freq,Er,LdB);     % Input impedance



% Plot the results on a smith chart (figure3 default)
smith(1,50);            % Plot Smith Chart at scale=1 and Zo=50 Ohms
smdrawc(Zload,50,'g-'); % Plot the load impedance Zload using Zo=50 Ohms
smdrawc(Z1a,50,'c-');   % Plot intermediate impedance Z1a using Z0=50 Ohms
smdrawc(Z2a,50,'m-');   % Plot intermediate impedance Z2a using Z0=50 Ohms
smdrawc(Zin,50,'r-');   % Plot the impedance Zin using Zo=50 Ohms 
smarker1(Zload,Freq,Zo,MidFreq,1);  % Put marker No.1 on Zload
smarker1(Z1a,Freq,Zo,MidFreq,2);    % Put marker No.2 on Z1a 
smarker1(Z2a,Freq,Zo,MidFreq,3);    % Put marker No.3 on Z2a
smarker1(Zin,Freq,Zo,MidFreq,4);    % Put marker No.4 on Zin
rlossc(Zin,Freq,Zo,'r');            % Plot return loss for Zin
