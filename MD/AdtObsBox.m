addpath(genpath('/home/nick/HDD/Work/Matlab-scripts'));

%%

clear all
DataDir='/afs/cern.ch/work/n/nbiancac/scratch0/2015_11_05/';

L=dir(DataDir);
L=L(3:end);
time_vec=[];
for ii=1:length(L)
    L(ii).time=datenum(L(ii).date,'dd-mm-yyyy HH:MM:SS');
    time_vec=[time_vec,L(ii).time];
end

[~,ind]=sort(time_vec);
L=L(ind);

last=L(end-1).name;
disp(last)
%%
close all
ResultDir='/home/nick/HDD/Work/CERN/MD/LHC/2015/4574/ADT/';
filename='MD754_ADTObsBox_20151105_102511_156821.h5';
%filename='MD754_ADTObsBox_20151105_035508_997396.h5';
flagshow='on';
flagsave=0;
for beam=[{'B2'}]
    for plane=[{'vertical'}];
       
        turn_vec=[1:1000:22001]; 
        delta=1000;

        B2H.data_new=[];
        for kk=1:length(turn_vec)-1
            turn1=turn_vec(kk);
            turn2=delta;
            try B2H.raw=h5read([DataDir,filename],['/',char(beam),'/',char(plane),],[1 turn1],[3564, turn2]); 
                disp([char(beam),' ',char(plane),' acq. turn ',num2str(turn1),'->',num2str(turn1+delta-1)]); 
                L=unwrap_ADT(B2H);
                B2H.data_new=[B2H.data_new,L.data];
                B2H.bucket=L.bucket;
            end
        end

        figure
        set(gcf,'visible',flagshow);
        stackedplot(B2H.data_new',4,5,4,'-b')
        ind=1:4:length(B2H.bucket);
        set(gca,'xtick',ind)
        set(gca,'xticklabel',B2H.bucket(ind))
        xlabel('bucket number')
        ylabel('Turns')
        zlabel('arb.units')
        title([char(beam),' ',char(plane)])
        
        if flagsave 
            name=[regexprep(filename,'.h5',''),'_',char(beam),'_',char(plane)];
            disp(name)
            s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='13';
            hgexport(gcf,'',s,'applystyle',true);
            hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');  
        end

        f_vec=[];fNAFF_vec=[];
        for ii=1:size(B2H.data_new,1)

            X=linspace(-0.5,0.5,length(B2H.data_new));
            freqabs=abs(fftshift(fft(B2H.data_new(ii,:))));
            ind=[X>0.25 & X< 0.35];

            freqabs_pos=freqabs(ind);
            freq_pos=X(ind);

            [~,ind_m]=max(freqabs_pos);
            f=freq_pos(ind_m);
            f_vec=[f_vec,f];
            
                    ind_rm=[X<0.295 & X> -0.295 | X>0.31 | X<-0.31];
                    spec=(fftshift(fft(B2H.data_new(ii,:))));
                    spec(ind_rm)=0;
                    newdata=real(ifft(fftshift(spec)));
                    
                    [avetune,errtune,tunes,phaseadv,amplit,numturns,frespec,nkick] = freqan(newdata',1,'all');
                    
            fNAFF_vec=[fNAFF_vec,avetune];
        end
        
        figure()
%         plot(B2H.bucket,f_vec-f_vec(1),'ok','markerfacecolor','k'); hold on;
        plot(B2H.bucket,fNAFF_vec,'or','markerfacecolor','r')
        title([char(beam),' ',char(plane)])
        ylabel('\Delta Q w.r.t. first train bunch')
        
        xlabel('Bucket number')
        grid on
       
        if flagsave 
            name=[regexprep(filename,'.h5',''),'_tune_train_',char(beam),'_',char(plane)];
            disp(name)
            s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='13';
            hgexport(gcf,'',s,'applystyle',true);
            hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');  
        end
    end
end
