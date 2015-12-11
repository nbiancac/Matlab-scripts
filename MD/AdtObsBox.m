addpath(genpath('/home/nick/HDD/Work/Matlab-scripts'));

%%
close all
clear all
% 4574 2x36b
DataDir='/afs/cern.ch/work/n/nbiancac/scratch0/2015_11_05/';
ResultDir='/home/nick/HDD/Work/CERN/MD/LHC/2015/4574/ADT/';
filename='ADTObsBox_20151105_015240_367084.h5'; %B1
% filename='ADTObsBox_20151105_015143_098178.h5'; %B2
QlimitsV=[0.294 0.3]; % for tune search in V
QlimitsH=[0.27 0.29]; % for tune search in H 

% 4577 72b
% DataDir='/afs/cern.ch/work/n/nbiancac/scratch0/2015_11_05/';
% ResultDir='/home/nick/HDD/Work/CERN/MD/LHC/2015/4577/ADT/';
% filename='ADTObsBox_20151105_092216_620029.h5'; %B1
% filename='ADTObsBox_20151105_092308_682746.h5'; %B2
% QlimitsV=[0.292 0.302]; % for tune search in V
% QlimitsH=[0.275 0.283]; % for tune search in H
nharm=50;
delta=200;
turn_vec=[1:delta:4001];
% turn_vec=[1:1000:22001]; 
% delta=1000;

flagshow='on';
flagsave=0;
flag_plot_fft=1;

ind_start=1; % start turn from injection
ind_end=450; % end turn FFT

 
        
for beam=[{'B1'},{'B2'}]
    for plane=[{'vertical'}]
       
        if strcmp(plane,'vertical'); Qlimits=QlimitsV; else Qlimits=QlimitsH; end
        ADT.data_new=[];
        for kk=1:length(turn_vec)-1
            turn1=turn_vec(kk);
            turn2=delta;
            try ADT.raw=h5read([DataDir,filename],['/',char(beam),'/',char(plane),],[1 turn1],[3564, turn2]); 
                disp([char(beam),' ',char(plane),' acq. turn ',num2str(turn1),'->',num2str(turn1+delta-1)]); 
                L=unwrap_ADT(ADT);
                if size(L.data,1)>size(ADT.data_new,1) % means we are over injecting
                    disp('Overinjection detected: reset matrix.')
                    ADT.data_new=[L.data];
                else
                    ADT.data_new=[ADT.data_new,L.data];
                end
                ADT.bucket=L.bucket;
            end
        end
        
        for kk=1:size(ADT.data_new,1)
            ADT.data_new(kk,:)=ADT.data_new(kk,:)-nanmean(ADT.data_new(kk,:));
        end
        
                
        
        
        if ~isempty(ADT.data_new)
            figure
            set(gcf,'visible',flagshow);
            stackedplot(ADT.data_new', 4 ,5,4,'-b')

            if length(ADT.bucket)>1
            ind=1:2:length(ADT.bucket);

            set(gca,'xtick',ind)
            set(gca,'xticklabel',ADT.bucket(ind))
            end
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
            disp('Getting index of injection')
            ADT.data_cut=[];
            ADT.bucket_cut=[];
            index_start_vec=[];
            for kk=1:size(ADT.data_new,1)
                index_start=find(diff(ADT.data_new(kk,:))>1e3);
                if ~isempty(index_start); index_start=index_start(1); 
                index_start_vec=[index_start_vec,index_start];
                end
            end
            
            index_start=min(index_start_vec);
            for kk=1:size(ADT.data_new,1)    
                ADT.data_cut=[ADT.data_cut;ADT.data_new(kk,index_start:end)];
                ADT.bucket_cut=[ADT.bucket_cut,ADT.bucket(kk)];
            end
            
            disp('remove pilot')
            ADT.data_cut(1,:)=[];
                ADT.bucket_cut(1)=[];

            figure;
            set(gcf,'visible',flagshow);
            stackedplot(ADT.data_cut', 4 ,5,4,'-b')

            if length(ADT.bucket_cut)>1
            ind=1:2:length(ADT.bucket_cut);

            set(gca,'xtick',ind);
            set(gca,'xticklabel',ADT.bucket_cut(ind));
            end
            xlabel('bucket number');
            ylabel('Turns');
            zlabel('arb.units');
            title([char(beam),' ',char(plane)]);

            if flagsave 
                name=[regexprep(filename,'.h5',''),'_',char(beam),'_',char(plane)];
                disp(name)
                s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='13';
                hgexport(gcf,'',s,'applystyle',true);
                hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');  
            end



            ADT.data_cut_fft=ADT.data_cut(:,ind_start:ind_end);
            
            
            figure;
            set(gcf,'visible',flagshow);
            stackedplot(ADT.data_cut_fft', 4 ,5,4,'-b');

            if length(ADT.bucket_cut)>1
                ind=1:2:length(ADT.bucket_cut);
                set(gca,'xtick',ind);
                set(gca,'xticklabel',ADT.bucket_cut(ind));
            end
            xlabel('bucket number');
            ylabel('Turns');
            zlabel('arb.units');
            title([char(beam),' ',char(plane)]);

            if flagsave 
                name=[regexprep(filename,'.h5',''),'_',char(beam),'_',char(plane)];
                disp(name)
                s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='13';
                hgexport(gcf,'',s,'applystyle',true);
                hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');  
            end

            %% damping time
           
%             figure(); semilogy(abs(hilbert(ADT.data_cut_fft(1,:)))./max(abs(hilbert(ADT.data_cut_fft(1,:)))),'ok')
% 
%             turn2fit=5:40;
%             data2fit=ADT.data_cut_fft(1,5:40);
%             p=polyfit(turn2fit,log(data2fit),1);
%             dampfit=p(1);
%             hold on
%             semilogy(turn2fit,exp(turn2fit.*dampfit),'r')
%             xlim([0 1500])
%             -1/dampfit
            
            
            %%
            
            f_vec=[];fNAFF_vec=[];
            for ii=1:size(ADT.data_cut_fft,1);
                
                ADT.data_cut_fft(ii,:)=ADT.data_cut_fft(ii,:)-mean(ADT.data_cut_fft(ii,:));
                X=linspace(-0.5,0.5,size(ADT.data_cut_fft,2));
                freqabs=abs(fftshift(fft(ADT.data_cut_fft(ii,:))));

                ind=[X>0. & X< 0.5];

                freqabs_pos=freqabs(ind);
                freq_pos=X(ind);

%                 figure(3)
%                 plot(X,freqabs);hold on;

                [~,ind_m]=max(freqabs_pos);
                f=freq_pos(ind_m);
                f_vec=[f_vec,f];


                
                newdata=ADT.data_cut_fft(ii,:);
%%
        
                [avetune,errtune,tunes,phaseadv,amplit,numturns,frespec,nkick] = freqan(newdata',nharm,'all');
                allfreq=abs(frespec(1:nharm,3));
                allamp=frespec(1:nharm,4);
                ind=[allfreq>Qlimits(1) &  allfreq<Qlimits(2)];
                freqlimits=allfreq(ind);
                amplimits=allamp(ind);
                [avetuneamp,indmax]=max(amplimits);
                avetune=freqlimits(indmax);
                if isempty(avetune); avetuneamp=nan; avetune=nan; end
                
                if flag_plot_fft==1
                    figure(121)
                    stem(allfreq,allamp,'.k'); hold on
                    set(gca,'yscale','log')
                    vline(Qlimits,'-r')
                    xlim([0.1 0.4])
                    plot(avetune, avetuneamp,'or'); hold off;
                end
                
                fNAFF_vec=[fNAFF_vec,avetune];
            end
            
            figure(2)
            plot(ADT.bucket_cut,fNAFF_vec,'or','markerfacecolor','r')
            title([char(beam),' ',char(plane)])
            ylim(Qlimits)
            if strcmp(char(plane),'horizontal');  ylabel(['Q_X']); else ylabel(['Q_Y']); end
            
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
%         pause
    end
end


%% Hilbert method

close all
clear all
% 4574 2x36b
DataDir='/afs/cern.ch/work/n/nbiancac/scratch0/2015_11_05/';
ResultDir='/home/nick/HDD/Work/CERN/MD/LHC/2015/4574/ADT/';
filename='ADTObsBox_20151105_015240_367084.h5'; %B1
% filename='ADTObsBox_20151105_015143_098178.h5';
QlimitsV=[0.294 0.3]; % for tune search in V
QlimitsH=[0.277 0.283]; % for tune search in H 

nharm=75;
delta=500;
turn_vec=[1:delta:4001];
% turn_vec=[1:1000:22001]; 
% delta=1000;

flagshow='on';
flagsave=0;
flag_plot_fft=1;

ind_start=1; % start turn from injection
ind_end=500; % end turn FFT

 
        
for beam=[{'B1'}]
    for plane=[{'horizontal'}]
       
        if strcmp(plane,'vertical'); Qlimits=QlimitsV; else Qlimits=QlimitsH; end
        ADT.data_new=[];
        for kk=1:length(turn_vec)-1
            turn1=turn_vec(kk);
            turn2=delta;
            try ADT.raw=h5read([DataDir,filename],['/',char(beam),'/',char(plane),],[1 turn1],[3564, turn2]); 
                disp([char(beam),' ',char(plane),' acq. turn ',num2str(turn1),'->',num2str(turn1+delta-1)]); 
                L=unwrap_ADT(ADT);
                if size(L.data,1)>size(ADT.data_new,1) % means we are over injecting
                    disp('Overinjection detected: reset matrix.')
                    ADT.data_new=[L.data];
                else
                    ADT.data_new=[ADT.data_new,L.data];
                end
                ADT.bucket=L.bucket;
            end
        end
        
        for kk=1:size(ADT.data_new,1)
            ADT.data_new(kk,:)=ADT.data_new(kk,:)-nanmean(ADT.data_new(kk,:));
        end
        
                
        
        
        if ~isempty(ADT.data_new)
            figure
            set(gcf,'visible',flagshow);
            stackedplot(ADT.data_new', 4 ,5,4,'-b')

            if length(ADT.bucket)>1
            ind=1:2:length(ADT.bucket);

            set(gca,'xtick',ind)
            set(gca,'xticklabel',ADT.bucket(ind))
            end
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
            disp('Getting index of injection')
            ADT.data_cut=[];
            ADT.bucket_cut=[];
            index_start_vec=[];
            for kk=1:size(ADT.data_new,1)
                index_start=find(diff(ADT.data_new(kk,:))>1e3);
                if ~isempty(index_start); index_start=index_start(1); 
                index_start_vec=[index_start_vec,index_start];
                end
            end
            
            index_start=min(index_start_vec);
            for kk=1:size(ADT.data_new,1)    
                ADT.data_cut=[ADT.data_cut;ADT.data_new(kk,index_start:end)];
                ADT.bucket_cut=[ADT.bucket_cut,ADT.bucket(kk)];
            end
            
            disp('remove pilot')
            ADT.data_cut(1,:)=[];
                ADT.bucket_cut(1)=[];

            figure;
            set(gcf,'visible',flagshow);
            stackedplot(ADT.data_cut', 4 ,5,4,'-b')
%%
close all
hmat=hilbert(ADT.data_cut(:,1:end));
 f_vec=[];frac_vec=[];
for ii=1:size(hmat,1)
            turns=1:size(hmat,2);
            phase=angle(hmat(ii,:)); 
            phase1=unwrap((phase),[],2)/2/pi;
            
            cv=diff(phase1)
            cv(cv<0)=cv(cv<0)+0.5
            phase2=cumsum(cv);
            
%             figure
%             plot(phase1); hold on;
%             plot(phase2,'r')
%             xlim([0 20])
            
            
            figure(1);
            plot(phase1,'-k'); hold on;
            plot(phase2,'r'); hold off
            xlim([0 20])
%             pause
%             gp=ginput(2)
            xlim1=5%;round(gp(1,1));
            xlim2=15%round(gp(2,1));
            
            
            xlim([0 20])

   
            turn2fit=xlim1:xlim2;
            phase2fit=phase2(turn2fit);

            figure(3)
            plot(turn2fit,phase2fit); hold on
            
            p=polyfit(turn2fit,(phase2fit),1);
            p(1)
            
            frac=p(1)-floor(p(1));
            
            if frac>=0.5
                frac=frac-0.5;
            end
            
            f_vec=[f_vec,frac];
            
            frac_vec=[frac_vec,frac];
end
figure(4)
plot(f_vec,'ok')

figure(5)
plot(frac_vec)
            
        end
    end
end
            