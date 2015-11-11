function CC=makeCCmodel(CrabDir,crabtype, plane, Ncavity, IP, spreadmode, FM,flag_rng, flagshow, betaCC, betav)
% function CC=makeCCmodel(CrabDir,crabtype, plane, Ncavity, IP, spreadmode, FM,flag_rng, flagshow, betaCC, betav)
%
% Function to collect all the HOM from Crab cavities starting from the
% table given by Rama & co. The scripts account for the crossing plane at
% the IPs: ie, since the crabs are designed with vertical crossing, in IP1
% we cross vertically and the crab is as in the table; in IP5 the cross is
% horizonthal and we need to rotate by 90deg the crab (and corresponding
% HOMs). The script account for beta, number of cavities, spread of the
% HOMs by 3MHz, and the possibility of having a feedback on the
% foundamental mode.

addpath(genpath('/afs/cern.ch/user/n/nbiancac/scratch0/Matlab-scripts/'));
%

HOM_table=dlmread([CrabDir,'HLLHC_Crab_cavities_',crabtype,'.txt'],'',1,0);    
disp(['Total ',num2str(Ncavity),' cavities @ beta=',num2str(betaCC),' in ',IP])

if plane~='z'
    if strcmp(IP,'IP1')
        cross_plane='y';
        disp([IP,' has ',cross_plane,' crossing: the CC is oriented in the same way.'])
        if plane==cross_plane
            mode_plane='y';
        else
            mode_plane='x';
        end

    elseif strcmp(IP,'IP5')
        cross_plane='x';
        disp([IP,' has ',cross_plane,' crossing: the CC is rotated.'])
        if plane==cross_plane
            mode_plane='y';
        else
            mode_plane='x';
        end
    end
else
    disp('Longitudinal modes: no differences at IPs.')
    mode_plane='z';
end


% initialize freq and impedance
freq_tot=10.^(1:1:10);      
Zt_tot=0*freq_tot;
   
    
if flag_rng==1; rng(1); end   % init the random generator         

disp('Calculating fundamental mode:')
i_mode=1; 
if strcmp(mode_plane,'z')
    R=HOM_table(i_mode,1);
    Q=HOM_table(i_mode,2);
    f=HOM_table(i_mode,3);
elseif strcmp(mode_plane,'x')
    R=HOM_table(i_mode,4)*betaCC/betav;
    Q=HOM_table(i_mode,5);
    f=HOM_table(i_mode,6);
elseif strcmp(mode_plane,'y')
    R=HOM_table(i_mode,7)*betaCC/betav;
    Q=HOM_table(i_mode,8);
    f=HOM_table(i_mode,9);
end
disp(['FM @ ',num2str(f/1e6),'MHz']);


if strcmp(FM,'RovQ') && ~strcmp(mode_plane,'z')
      disp('FM damped with R/Q=1.')
      R=R/Q;
      Q=1;
      disp('No spread assumed on FM.');
      f_spread=zeros(1,8); % uniform spread +/- 3MHz
      for ii=1:Ncavity; % sum up the FM

        f=f+f_spread(ii);
        N=20;
        D=f/2/Q;
        freq_mode=f;
        for nn=10.^((0:2:11)/2)
                frmin=min(f);frmax=max(f);
                freq_mode=[freq_mode,frmin-nn*D:(nn*D)/N:frmin];
                freq_mode=[freq_mode,frmax:(nn*D)/N:frmax+nn*D];
        end

        freq_mode=sort(unique(freq_mode(freq_mode>0)));
        if strcmp(mode_plane,'z')
            Zt_mode=long_resonator(freq_mode,Q,f,R);
        else
            Zt_mode=transverse_resonator(freq_mode,Q,f,R);
        end
        freq_tot_new=[freq_tot,freq_mode];
        freq_tot_new=unique(sort(freq_tot_new(freq_tot_new>0)));
        Zt_tot_new=interp1(freq_tot,Zt_tot,freq_tot_new);
        Zt_mode=interp1(freq_mode,Zt_mode,freq_tot_new);
        Zt_tot=Zt_tot_new+Zt_mode;
        freq_tot=freq_tot_new;

      end
elseif strcmp(FM,'Feedback') && ~strcmp(mode_plane,'z')
      disp('FM with a feedback on it.')
      L=dlmread([CrabDir,'1stCCmodewithfeedback_clean2.dat']);
      [freqFM_int,ReFM_int,ImFM_int]=deal(L(:,1),L(:,2),L(:,3));
      freq_tot=unique(sort([freq_tot,(freqFM_int')]));
      Zt_tot=interp1(freqFM_int,ReFM_int+1i*ImFM_int,freq_tot)*Ncavity;
elseif strcmp(FM,'None') || strcmp(mode_plane,'z')
      disp('no FM.')
      Zt_tot=0*freq_tot;
end 


disp(['Calculating HOMs for ',IP,' in ',plane,'-plane']);
disp(['that corresponds to the HOM in the ',mode_plane,'-plane of the crab cavity.']);
freq_base=sort([10.^(1:1/1:10),1e6:1e6:5e9]);      
for i_mode=2:size(HOM_table,1)

    
    if strcmp(mode_plane,'z')
        R=HOM_table(i_mode,1);
        Q=HOM_table(i_mode,2);
        f=HOM_table(i_mode,3);
    elseif strcmp(mode_plane,'x')
        R=HOM_table(i_mode,4)*betaCC/betav;
        Q=HOM_table(i_mode,5);
        f=HOM_table(i_mode,6);
    elseif strcmp(mode_plane,'y')
        R=HOM_table(i_mode,7)*betaCC/betav;
        Q=HOM_table(i_mode,8);
        f=HOM_table(i_mode,9);
    end
    disp(['HOM: ', num2str(i_mode),'/',num2str(size(HOM_table,1)),' @ ',num2str(f/1e6),'MHz'])
    
    if strcmp(spreadmode,'On')
        f_spread=-3e6+6e6*rand(1,8); % uniform spread +/- 3MHz
    else
        f_spread=0*rand(1,8); % uniform spread +/- 3MHz
    end
    for ii=1:Ncavity; % sum up the HOM for each crab on the corrsponding plane
        f=f+f_spread(ii);
        N=20; % hardcoded: points per resonance decade
        D=f/2/Q;
        freq_mode=f;
        for nn=10.^((0:2:11)/2)
                frmin=min(f);frmax=max(f);
                freq_mode=[freq_mode,frmin-nn*D:(nn*D)/N:frmin];
                freq_mode=[freq_mode,frmax:(nn*D)/N:frmax+nn*D];
        end
        freq_mode=sort(unique([freq_base,freq_mode(freq_mode>0)]));
        if strcmp(mode_plane,'z')
            Zt_mode=long_resonator(freq_mode,Q,f,R);
        else
            Zt_mode=transverse_resonator(freq_mode,Q,f,R);
        end
        freq_tot_new=[freq_tot,freq_mode];
        freq_tot_new=unique(sort(freq_tot_new(freq_tot_new>0)));
        Zt_tot_new=interp1(freq_tot,Zt_tot,freq_tot_new);
        Zt_mode=interp1(freq_mode,Zt_mode,freq_tot_new);
        Zt_tot=Zt_tot_new+Zt_mode;
        freq_tot=freq_tot_new;
    end

end

ind=isnan(Zt_tot);
Zt_tot(ind)=[];
freq_tot(ind)=[];

if strcmp(flagshow,'on')
    figure(2);
    set(gcf,'visible',flagshow);
    semilogy(freq_tot,real(Zt_tot),'-'); hold on;
    xlim([0 2e9])
end

CC.Z=Zt_tot;
CC.freq=freq_tot;
end

