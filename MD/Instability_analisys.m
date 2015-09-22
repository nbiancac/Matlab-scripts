% addpath(genpath('/afs/cern.ch/user/n/nbiancac/scratch0/Matlab-scripts/'));

addpath(genpath('/home/nick/HDD/Work/Matlab-scripts'));

% mainDir='/home/nick/HDD/Dropbox/CERN/MD/4254/';
fill='29082015-ACdipole';
mainDir=['/home/nick/HDD/Work/CERN/MD/LHC/2015/',fill,'/'];
DataDir=[mainDir,'Data/'];
ResultDir=[mainDir,'Result/'];
CollDir=[DataDir,'Collimators/'];

% machine parameters
[~,~,~,E0]=particle_param('proton');
machine=LHC_param(E0,450e9,'Nominal LHC');
disp([char(machine.scenario),' @ ',machine.Estr])
% beam and plane
beam='2';
plane='V';
comment=''; % for BBQ

flagsave=0;

disp(['Analysis for B',beam,plane,' on fill ',fill])
%% Get BBQ raw

typeBBQ='HS';
if strcmp(typeBBQ,'HS')
    variable=['LHC.BQBBQ.CONTINUOUS_HS.B',beam,'_ACQ_DATA_',plane,comment];
else
    variable=['LHC.BQBBQ.CONTINUOUS.B',beam,'_ACQ_DATA_',plane,comment];
end
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');

delta_time_BBQ=time(2)-time(1);
turns_interval=round(eval(datestr(delta_time_BBQ,'SS.FFF'))*machine.f0); % timing between two acquisitions
turns_acq=size(L,2); % number of acquired turns
turns_gap=turns_interval-turns_acq; % >0 gap, <0 overlap
if turns_gap>0
    L=[L,zeros(size(L,1),turns_gap)];
else
    L(:,end+turns_gap:end)=[];
end

BBQ={};
BBQ.data=reshape(transp(L),1,[]); % after cleaning overlaps or gap we connect all the turns.
BBQ.turns=1:length(BBQ.data);
BBQ.time=linspace(time(1),time(end),length(BBQ.turns));

figure(1)
plot(BBQ.time,BBQ.data)
N=10;
xtime=linspace(BBQ.time(1),BBQ.time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'));
set(gca,'xtick',(xtime));
rotateticklabel(gca,90,10);
%% Get variables
% BBQG & BBQHS

BBQG={};
if strcmp(plane,'H')
variable=['LHC.BQBBQ.CONTINUOUS.B',beam,'_EIGEN_AMPL_1']; else
variable=['LHC.BQBBQ.CONTINUOUS.B',beam,'_EIGEN_AMPL_2']; 
end
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
BBQG.data=L;
BBQG.time=time;
turns_interval=round(eval(datestr(time(2)-time(1),'SS.FFF'))*machine.f0); % timing between two acquisitions
BBQG.turns=turns_interval*[1:length(time)];

BBQHS={};
if strcmp(plane,'H')
variable=['LHC.BQBBQ.CONTINUOUS_HS.B',beam,'_EIGEN_AMPL_1']; else
variable=['LHC.BQBBQ.CONTINUOUS_HS.B',beam,'_EIGEN_AMPL_2']; 
end
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
BBQHS.data=L;
BBQHS.time=time;
turns_interval=round(eval(datestr(time(2)-time(1),'SS.FFF'))*machine.f0); % timing between two acquisitions
BBQHS.turns=turns_interval*[1:length(time)];

% FBCT

FBCT={};
variable=['LHC.BCTFR.A6R4.B',beam,'_BUNCH_INTENSITY'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
FBCT.dataall=L;
for ind=1:size(FBCT.dataall,2)
    FBCT.dataall_rel(:,ind)=FBCT.dataall(:,ind)./max(FBCT.dataall(:,ind));
end
turns_interval=round(eval(datestr(time(2)-time(1),'SS.FFF'))*machine.f0); % timing between two acquisitions
FBCT.turns=turns_interval*[1:length(time)];
FBCT.time=time;

% BCT
BCT={};
variable=['LHC.BCTDC.A6R4.B',beam,'_BEAM_INTENSITY'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
BCT.data=L;
turns_interval=round(eval(datestr(time(2)-time(1),'SS.FFF'))*machine.f0); % timing between two acquisitions
BCT.turns=turns_interval*[1:length(time)];
BCT.time=time;

% ROF
ROF={};
variable=['RPMBB.RR13.ROF.A81B',beam,'_I_MEAS'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
ROF.data=L;
turns_interval=round(eval(datestr(time(2)-time(1),'SS.FFF'))*machine.f0); % timing between two acquisitions
ROF.turns=turns_interval*[1:length(time)];
ROF.time=time;

% BETAs

BETA1={};BETA5={};BETA8={};
variable=['TCDD.4L2_EXPACQ_BETASTAR_IP1'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
BETA1.data=L;
turns_interval=round(eval(datestr(time(2)-time(1),'SS.FFF'))*machine.f0); % timing between two acquisitions
BETA1.turns=turns_interval*[1:length(time)];
BETA1.time=time;

variable=['TCDD.4L2_EXPACQ_BETASTAR_IP5'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
BETA5.data=L;
turns_interval=round(eval(datestr(time(2)-time(1),'SS.FFF'))*machine.f0); % timing between two acquisitions
BETA5.turns=turns_interval*[1:length(time)];
BETA5.time=time;

variable=['TCDD.4L2_EXPACQ_BETASTAR_IP8'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
BETA8.data=L;
turns_interval=round(eval(datestr(time(2)-time(1),'SS.FFF'))*machine.f0); % timing between two acquisitions
BETA8.turns=turns_interval*[1:length(time)];
BETA8.time=time;

% BSRT
BSRT={};
if strcmp(beam,'2'); section='L'; else section='R'; end
variable=['LHC.BSRT.5',section,'4.B',beam,'_GATE_DELAY'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
turns_interval=round(eval(datestr(time(2)-time(1),'SS.FFF'))*machine.f0); % timing between two acquisitions
time_interval=time(2)-time(1);
BSRT.bucket=L;

BSRT.turns=[]; BSRT.time=[];

for ii=1:length(time)
    BSRT.turns=[BSRT.turns,(ii-1)*turns_interval+[linspace(0,turns_interval,size(BSRT.bucket,2)+1)]];
    BSRT.time=[BSRT.time,time(ii)+[0:(size(BSRT.bucket,2)-1)]*(time_interval/size(BSRT.bucket,2))];

end

BSRT.bucket=1+[reshape(BSRT.bucket',1,[])];


variable=['LHC.BSRT.5',section,'4.B',beam,'_FIT_SIGMA_H'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
BSRT.sigmaH=reshape(L',[],1); 


variable=['LHC.BSRT.5',section,'4.B',beam,'_FIT_SIGMA_V'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
BSRT.sigmaV=reshape(L',[],1); 
 
variable=['LHC.BSRT.5',section,'4.B',beam,'_LSF_H'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
BSRT.LSFH=L; 

variable=['LHC.BSRT.5',section,'4.B',beam,'_LSF_V'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
BSRT.LSFV=L;

variable=['LHC.BSRT.5',section,'4.B',beam,'_BETA_H'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
BSRT.BETAH=L; 

variable=['LHC.BSRT.5',section,'4.B',beam,'_BETA_V'];
[time,L]=read_timber_data([DataDir,variable,'.dat'],'yyyy-mm-dd HH:MM:SS.FFF');
BSRT.BETAV=L;

BSRT.eH=((BSRT.sigmaH).^2-BSRT.LSFH^2)/(BSRT.BETAH)*(machine.beta*machine.gamma);
BSRT.eV=((BSRT.sigmaV).^2-BSRT.LSFV^2)/(BSRT.BETAV)*(machine.beta*machine.gamma);


%% Get collimators
comment_scenario='Scan_TDI_closed/';
L=dir([CollDir,comment_scenario,'*B',beam,'*LVDT_GU.dat']);
%%
Coll=[];

for ii=1:length(L)
    Coll(ii).name={};
    Coll(ii).data={};
    disp([L(ii).name]);
    namecoll=regexprep(L(ii).name,[{'_MEAS_LVDT_GU.dat'},{'TCTP'},{'TCSP.A4R6'},{'TCTV'}],[{''},{'TCT'},{'TCSG.4R6'},{'TCTVA'}]);
    [time,data]=read_timber_data([CollDir,comment_scenario,L(ii).name],'yyyy-mm-dd HH:MM:SS.FFF');
    Coll(ii).name=namecoll;
    Coll(ii).data=data;
    Coll(ii).time=time;
end
%%
if strcmp(beam,'1')
    L=dir([DataDir,'TCLIA.4R2_MEAS_LVDT_GU.dat']);
    L=[L,dir([DataDir,'TDI.4L2_MEAS_LVDT_GU.dat'])];
    mm=length(Coll);
    for ii=1:length(L)
        disp([L(ii).name]);
        namecoll=regexprep(L(ii).name,[{'_MEAS_LVDT_GU.dat'},{'TCTP'},{'TCSP.A4R6'},{'TCTV'}],[{''},{'TCT'},{'TCSG.4R6'},{'TCTVA'}]);
        [time,data]=read_timber_data([CollDir,comment_scenario,L(ii).name],'yyyy-mm-dd HH:MM:SS.FFF');
        Coll(mm+ii).name=namecoll;
        Coll(mm+ii).data=data;
        Coll(mm+ii).time=time;
    end
else
    L=dir([CollDir,comment_scenario,'TCLIA.4L8_MEAS_LVDT_GU.dat']);
    L=[L,dir([CollDir,comment_scenario,'TDI.4R8_MEAS_LVDT_GU.dat'])];
    mm=length(Coll);
    for ii=1:length(L)
        disp([L(ii).name]);
        namecoll=regexprep(L(ii).name,[{'_MEAS_LVDT_GU.dat'},{'TCTP'},{'TCSP.A4R6'},{'TCTV'}],[{''},{'TCT'},{'TCSG.4R6'},{'TCTVA'}]);
        [time,data]=read_timber_data([CollDir,comment_scenario,L(ii).name],'yyyy-mm-dd HH:MM:SS.FFF');
        Coll(mm+ii).name=namecoll;
        Coll(mm+ii).data=data;
        Coll(mm+ii).time=time;
    end
end


%%
CollDataDir=['/afs/cern.ch/user/n/nbiancac/ln_delphi/PYTHON_codes_and_scripts/LHC_impedance_and_scripts/Coll_settings/'];
file1='LHC_ft_6.5TeV_B1.txt';
disp(['Making table with ',file1,' model.'])

[names1,gaps1,betx1,bety1,sigma1]=read_coll_file([CollDataDir,file1],'all','nounique');

count=0;
for ii=1:length(names1)
    disp(char(names1(ii)))
    for jj=1:length(Coll)
        if strcmp(char(Coll(jj).name),char(names1(ii)))
            count=count+1;
            disp([char(Coll(jj).name),' = ',char(names1(ii))])
        elseif (strncmp(char(Coll(jj).name),'TDI',3)) && (strncmp(char(names1(ii)),'TDI',3))
            count=count+1;
            disp([char(Coll(jj).name),' = ',char(names1(ii))])
        elseif (strncmp(char(Coll(jj).name),'TCDQ',4)) && (strncmp(char(names1(ii)),'TCDQ',4))
            count=count+1;
            disp([char(Coll(jj).name),' = ',char(names1(ii))])
        end
    end
end
%%
disp(['Collimators in table: ',num2str(length(names1))]);
disp(['Collimators from Timber: ',num2str(length(Coll))]);
disp(['Collimators tabbed from Timber to table: ',num2str(count)]);


%%
time_sample=Coll(1).time(1);
hg_vec=[];
for ii=1:length(Coll)
    hg_vec=[hg_vec,interp1(Coll(ii).time,Coll(ii).data,time_sample)/2];
end
%%
swap_coll_file([CollDataDir,file1],[CollDataDir,'prova'],'Halfgap[m]', hg_vec)

%% FBCT
close all
flagsave=0;
flagwrite=1;
flagshow='on';


filled_buckets=find(FBCT.dataall(floor(end/2),:)>5e10);
pilot_buckets=find(FBCT.dataall(floor(end/2),:)>1e9 & FBCT.dataall(floor(end/2),:)<5e10);

Mbunches=numel(filled_buckets);
disp(['M=',num2str(Mbunches)])
index=[filled_buckets];
FBCT.dataall(FBCT.dataall==0)=nan;
maxfbct=nanmax(FBCT.dataall,[],1)/1e11;

% Filling scheme
fill_pattern=diff(find((maxfbct>.1)));
fill_pattern(end+1)=1e6;
sum=0;
fill_scheme=[];
for ii=1:length(fill_pattern)
    if fill_pattern(ii)==1
        sum=sum+1;
    else
        fill_scheme=[fill_scheme,sum];
        sum=0;
    end
end

ind=find(fill_scheme>1);
fill_scheme(ind)=fill_scheme(ind)+1;
disp(fill_scheme);
dlmwrite([ResultDir,'filling_scheme.txt'],fill_scheme);


% All FBCT and buckets
f=figure(13);
set(gcf,'visible',flagshow)
h=imagesc(1:size(FBCT.dataall,2),FBCT.time,FBCT.dataall_rel);
caxis([0.9 1])
xlim([0 3565])
ylabel('Time'); xlabel('# bucket');
h=colorbar;
title(['Fill:',fill,' FBCT B',beam])
N=10;
ytime=linspace(FBCT.time(1),FBCT.time(end),N);
set(gca,'yticklabel',datestr(ytime,'HH:MM'))
set(gca,'ytick',(ytime));
set(gcf,'position', [407   268   794   384])

if flagsave
    name=['FBCT_all_B',beam,plane];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='8';
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end
%

% FBCT for chosen index
FBCT_sel={};
FBCT_sel.turns=FBCT.turns;
FBCT_sel.time=FBCT.time;
FBCT_sel.data=[];
FBCT_sel.index=index;
col_vec=distinguishable_colors(length(index));
leg_vec=[];
for ii=1:length(index)
    figure(15)
    set(gcf,'visible',flagshow)
    plot(FBCT.time,FBCT.dataall(:,index(ii)),'color',col_vec(ii,:)); hold on;
    leg_vec=[leg_vec,{['# ',num2str(index(ii))]}];
    FBCT_sel.data=[FBCT_sel.data,FBCT.dataall(:,index(ii))];
    if flagwrite==1
        dlmwrite([DataDir,'intensity_b',num2str(index(ii)),'.dat'],[FBCT.time',FBCT.dataall(:,index(ii))],'delimiter', '\t', 'precision', 16);
    end
end
hold off;
FBCT_sel.leg=leg_vec;
if Mbunches<10
    legend(leg_vec,'location','southwest');
else
    legend('All filled buckets','location','southwest');
end
ylabel('N_b [ppb]');
title(['Fill:',fill,' FBCT B',beam])
grid on
N=10;
xtime=linspace(FBCT.time(1),FBCT.time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'))
set(gca,'xtick',(xtime))
rotateticklabel(gca,90,6)
set(gcf,'position', [407   268   794   384])

if flagsave
    name=['FBCT_vs_selected_bucket_B',beam,plane];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='8';
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end

% relative FBCT for chosen index
for ii=1:length(index)
    figure(17)
    set(gcf,'visible',flagshow)
    plot(FBCT.time,FBCT.dataall(:,index(ii))/FBCT.dataall(1,index(ii)),'color',col_vec(ii,:)); hold on;
end
hold off;
if Mbunches<10
    legend(leg_vec,'location','southwest');
else
    legend('All filled buckets','location','southwest');
end
ylabel('N_b/N^0_b');
title(['Fill:',fill,' FBCT B',beam,' M=',num2str(Mbunches)])
grid on
N=10;
xtime=linspace(FBCT.time(1),FBCT.time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'))
set(gca,'xtick',(xtime))
rotateticklabel(gca,90,6)
set(gcf,'position', [407   268   794   384])

if flagsave
    name=['FBCTperc_vs_selected_bucket_B',beam,plane];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='8';
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end

% min/max FBCT
figure(16)
set(gcf,'visible',flagshow)
stem(1:size(FBCT.dataall,2),max(FBCT.dataall,[],1)/1e11,'.r'); hold on;
stem(1:size(FBCT.dataall,2),min(FBCT.dataall,[],1)/1e11,'.k'); hold off
legend('before (max N_b)','after (min N_b)','location','southwest')
title(['Fill:',fill,' FBCT B',beam,' M=',num2str(Mbunches)])
xlim([1,3564])
xlabel('# bucket'); ylabel('max N_b [10^{11} ppb]');
set(gcf,'position', [407   268   794   384])

if flagsave
    name=['FBCT_vs_buckets_B',beam,plane];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='8';
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end


%% FBCT & ROF
flagsave=1;
close all
figure(1)
h=plot(FBCT_sel.time,FBCT_sel.data); 
h2=axisyy(ROF.time,ROF.data,'-m');
h3=axisyy(BBQG.time,BBQG.data,'--','color',[.5 .5 .5]);
h4=axisyy(BBQHS.time,BBQHS.data,'-g');
h2.YLabel='ROF';
h3.YLabel='BBQ gated';
h4.YLabel='BBQ HS';
legend(h,num2str(filled_buckets'),'location','southwest')
title(['B',beam,plane])
N=10;
xtime=linspace(BBQHS.time(1),BBQHS.time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'))
set(gca,'xtick',(xtime))
rotateticklabel(gca,90,10)
set(gcf,'position', [407   268   794   384])

if flagsave
    name=['FBCT_vs_BBQ_B',beam,plane];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='8';
    hgexport(gcf,'',s,'applystyle',true);
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end


%% FBCT & BBQ
flagsave=1;

nevery=1;

figure(1)
if Mbunches<10
    h1=plot(FBCT_sel.time,FBCT_sel.data./1e11); hold on
else
    h1=plot(FBCT_sel.time,sum(FBCT_sel.data,2)./1e11); hold on
end
h2=axisyy(BBQHS.time(1:nevery:end),BBQHS.data(1:nevery:end));

ylabel(['FBCT B',beam,plane,' [1e11 ppb]']);
h2.YLabel=['BBQ HS data B',beam,plane];

legend(leg_vec,'location','southwest');
N=10;
xtime=linspace(BBQHS.time(1),BBQHS.time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'))
set(gca,'xtick',(xtime))
rotateticklabel(gca,90,10)
set(gcf,'position', [407   268   794   384])

if flagsave
    name=['FBCT_B',beam,plane];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='8';
    hgexport(gcf,'',s,'applystyle',true);
%     saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end

%% BCT & BBQ
flagsave=0;

[BBQ,BCT]=interp_timber_data(BBQ,BCT,1); 

nevery=100;
figure(2)
[h]=plotyy(BBQ.turns(1:nevery:end),BBQ.data(1:nevery:end),BCT.turns,BCT.data);
ylabel(h(1),['BBQ data B',beam,plane])
ylabel(h(2),['BCT B',beam,plane])
xlabel('Turns')
set(gcf,'position', [407   268   794   384])

if flagsave
    name=['BCT_B',beam,plane,'_',comment];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='15';
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end

%% ROF & BBQ

figure(3)
[h]=plot(BBQ.time,BBQ.data);
h2=axisyy(ROF.time,ROF.data,'-k');
xlabel('Turns')
ylabel('BBQ [arb.units]')
h2.YLabel=('ROF [A]');
set(gcf,'position', [407   268   794   384])

title(['Fill:',fill,' B',beam,plane,'_',comment])
if flagsave
    name=['ROF_vs_BBQ_B',beam,plane];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='15';
    hgexport(gcf,'',s,'applystyle',true);
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end

%% BCT & ROF

figure(3)
plot(ROF.time,ROF.data);
xlabel('time')
ylabel('ROF [A]')
h2=axisyy(BCT.time,BCT.data,'-r');
h2.YLabel='BCT [p+]';


%% BETASTAR & BBQ
flagsave=1;
figure(3)
[h]=plot(BBQHS.time,BBQHS.data);
addaxis(BETA1.time,BETA1.data);
addaxis(BETA5.time,BETA5.data);
addaxis(BETA8.time,BETA8.data);
addaxislabel(1,'BBQ [arb.units]')
addaxislabel(2,'IP1 \beta^* [m]')
addaxislabel(3,'IP5 \beta^* [m]')
addaxislabel(4,'IP8 \beta^* [m]')
title(['Fill:',fill,' B',beam,plane])

N=10;
xtime=linspace(BBQHS.time(1),BBQHS.time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'))
set(gca,'xtick',(xtime))
rotateticklabel(gca,90,8);
xlim([BBQHS.time(1),BBQHS.time(end)]);
set(gcf,'position', [407   268   794   384])

if flagsave
    name=['BETASTAR_B',beam,plane];
    s=hgexport('readstyle','PRSTAB-6pt'); s.FixedFontSize='15';
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end

%% BSRT
close all

flagsave=0;
flagshow='on';

index=filled_buckets;

figure
set(gcf,'visible',flagshow);
h1=plot(BSRT.time,BSRT.bucket,'.-');
h2=axisyy(BSRT.time,BSRT.eV,'.-k');
h3=axisyy(BSRT.time,BSRT.eH,'.-m');
h4=axisyy(BBQHS.time,BBQHS.data,'-g');

ylabel('# bucket')
h2.YLabel='\epsilon_V [um]';
h3.YLabel='\epsilon_H [um]';
h4.YLabel='BBQ-HS [a.u.]';

N=10;
xtime=linspace(BBQHS.time(1),BBQHS.time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'));
set(gca,'xtick',(xtime));
rotateticklabel(gca,90,10);
set(gcf,'position', [407   268   794   384])
title(['Fill:',fill,' BBQHS  B',beam,plane])

if flagsave
    comment='';
    name=['BSRT_',comment,'B',beam,plane];
    s=hgexport('readstyle','PRSTAB-6pt'); s.FixedFontSize='8';
    hgexport(gcf,'',s,'applystyle',true);
%     saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end

%% seleted buckets

ind=[filled_buckets];

BSRT_sel={};
BSRT_sel.turns=[];
BSRT_sel.time=[];
BSRT_sel.eV=[];BSRT_sel.eH=[];BSRT_sel.bucket=[];
for ii=ind
   
    ind_eV=find(BSRT.bucket==ii);
    if ~isempty(ind_eV)
        BSRT_sel.turns=BSRT.turns;
        BSRT_sel.time=(BSRT.time);
        databsrt=interp1(BSRT.time(ind_eV),BSRT.eV(ind_eV),(BSRT.time),'linear');
        BSRT_sel.eV=[BSRT_sel.eV,databsrt'];
        BSRT_sel.bucket=[BSRT_sel.bucket,ii];
        
        ind_eH=find(BSRT.bucket==ii);
        databsrt=interp1(BSRT.time(ind_eH),BSRT.eH(ind_eH),(BSRT.time),'linear');
        BSRT_sel.eH=[BSRT_sel.eH,databsrt'];
    else
        disp(['Bucket ',num2str(ii),' not scanned.'])
    end
end


% BSRT vs bucket vs Squeeze vs BBQ
col_vec=distinguishable_colors(length(ind));
h1_vec=[]; leg1_vec=[];
h2_vec=[]; leg2_vec=[];
figure
set(gcf,'visible',flagshow);

for ii=1:length(ind)
        h1=plot(BSRT_sel.time,BSRT_sel.eV(:,ii),':','color',col_vec(ii,:)); hold on;
        h2=plot(BSRT_sel.time,BSRT_sel.eH(:,ii),'-','color',col_vec(ii,:)); hold on;
        h1_vec=[h1_vec,h1];
        h2_vec=[h2_vec,h2];
        leg1_vec=[leg1_vec,{['\epsilon_V #',num2str(BSRT_sel.bucket(ii))]}];
        leg2_vec=[leg2_vec,{['\epsilon_H #',num2str(BSRT_sel.bucket(ii))]}];
end

h2=axisyy(BETA1.time,BETA1.data);
h3=axisyy(BBQHS.time,BBQHS.data,'-');
if Mbunches<10
    legend([h1_vec,h2_vec],[leg1_vec,leg2_vec],'location','northwest');
else
    legend('All filled buckets','location','northwest');
end
N=10;
xtime=linspace(BBQHS.time(1),BBQHS.time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'))
set(gca,'xtick',(xtime))
rotateticklabel(gca,90,10);
set(gcf,'position', [407   268   794   384])

ylabel('\epsilon [um]')
h2.YLabel=('\beta* IP1 [m]');
h3.YLabel=('BBQ HS [a.u.]');

title(['Fill:',fill,' BBQHS  B',beam,plane])

if flagsave
    comment='';
    name=['BSRT_selected_',comment,'B',beam,plane];
    s=hgexport('readstyle','PRSTAB-6pt'); s.FixedFontSize='8';
    hgexport(gcf,'',s,'applystyle',true);
%     saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end

% BSRT BBQHS BBQG INT OCT

col_vec=distinguishable_colors(length(ind));
h1_vec=[]; leg1_vec=[];
h2_vec=[]; leg2_vec=[];
figure
set(gcf,'visible',flagshow);
smooth_wind=350;
for ii=1:length(ind)
        dataV=BSRT_sel.eV(:,ii);
        timeV=BSRT_sel.time(~isnan(dataV));
        dataV=dataV(~isnan(dataV));
        dataH=BSRT_sel.eH(:,ii);
        timeH=BSRT_sel.time(~isnan(dataH));
        dataH=dataH(~isnan(dataH));
        h1=plot(timeV,fastsmooth(dataV,smooth_wind),':','color',col_vec(ii,:)); hold on;
        h2=plot(timeH,fastsmooth(dataH,smooth_wind),'-','color',col_vec(ii,:)); hold on;
        h1_vec=[h1_vec,h1];
        h2_vec=[h2_vec,h2];
        leg1_vec=[leg1_vec,{['\epsilon_V #',num2str(BSRT_sel.bucket(ii))]}];
        leg2_vec=[leg2_vec,{['\epsilon_H #',num2str(BSRT_sel.bucket(ii))]}];
end

h3=axisyy(BBQG.time,BBQG.data,'YLim',[0,max(BBQHS.data)],'-','color',[.5 .5 .5]);
h4=axisyy(BBQHS.time,BBQHS.data,'-k');
h5=axisyy(ROF.time,ROF.data,'-m');

if Mbunches<10
for ii=1:length(ind)
    if ii==(1)
        h6=axisyy(FBCT_sel.time,FBCT_sel.data(:,ii)/1e11,'YLim',[0 1.5],'-','color',col_vec(ii,:)); hold on;
    else
        h6.plot(FBCT_sel.time,FBCT_sel.data(:,ii)/1e11,'-','color',col_vec(ii,:));
    end
end
else
    clear sum
    h6=plot(FBCT_sel.time,sum(FBCT_sel.data,2)/1e11);
end
if Mbunches<10
    legend([h1_vec,h2_vec],[leg1_vec,leg2_vec],'location','northwest');
else
    legend('All filled buckets','location','northwest');
end

N=20;
xtime=linspace(BBQHS.time(1),BBQHS.time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'))
set(gca,'xtick',(xtime))
rotateticklabel(gca,90,8);
title(['Fill:',fill,' B',beam,plane])
set(gcf,'position', [   387   145   912   507])

ylabel('\epsilon [um]')
h3.YLabel='BBQ G [arb. units]';
h4.YLabel='BBQ HS [arb. units]';
h5.YLabel='ROF [A]';
h6.YLabel='FBCT [10^{11} ppb]';


if flagsave
    comment='';
    name=['BSRT_FBCT_selected_',comment,'B',beam,plane];
    s=hgexport('readstyle','PRSTAB-6pt'); s.FixedFontSize='7';
    hgexport(gcf,'',s,'applystyle',true);
%     saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end

%%
flagsave=1;

figure()
plot(BBQHS.time,BBQHS.data,'-k'); hold on;
axisyy(BSRT.time,BSRT.bucket,'-','color',[.5 .5 .5]); hold on;
axisyy(BBQHS.time,BBQHS.data,'-k');

N=10;
xtime=linspace(BBQHS.time(1),BBQHS.time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'))
set(gca,'xtick',(xtime))
rotateticklabel(gca,90,10);

p=ginput(2);
xlim1=p(1,1);
xlim2=p(2,1);


figure
set(gcf,'visible',flagshow);
interval=find(BSRT_sel.time>xlim1 & BSRT_sel.time<xlim2);
for ii=1:length(ind)
        h1=stem(ind(ii),(BSRT_sel.eV(interval(1),ii)),'.-r'); hold on;
        h2=stem(ind(ii),(BSRT_sel.eV(interval(end),ii)),'.k'); hold on;
end

legend([h1,h2],datestr(xlim1,'HH:MM:SS'),datestr(xlim2,'HH:MM:SS'),'location','northwest'); legend boxoff


title(['Fill:',fill,' BSRT B',beam,' M=',num2str(Mbunches)])
xlim([1,3564])
xlabel('# bucket'); ylabel('\epsilon_V [um]');

if flagsave
    name=['BSRT_vs_buckets_zoom_B',beam,'V'];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='8';
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end

figure
set(gcf,'visible',flagshow);
interval=find(BSRT_sel.time>xlim1 & BSRT_sel.time<xlim2);
for ii=1:length(ind)
        h1=stem(ind(ii),(BSRT_sel.eH(interval(1),ii)),'.-r'); hold on;
        h2=stem(ind(ii),(BSRT_sel.eH(interval(end),ii)),'.k'); hold on;
end

legend([h1,h2],datestr(xlim1,'HH:MM:SS'),datestr(xlim2,'HH:MM:SS'),'location','northwest'); legend boxoff


title(['Fill:',fill,' BSRT B',beam,' M=',num2str(Mbunches)])
xlim([1,3564])
xlabel('# bucket'); ylabel('\epsilon_H [um]');

if flagsave
    name=['BSRT_vs_buckets_zoom_B',beam,'H'];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='8';
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end


figure
set(gcf,'visible',flagshow);
interval=find(FBCT_sel.time>xlim1 & FBCT_sel.time<xlim2);
for ii=1:length(ind)
        h1=stem(ind(ii),(FBCT_sel.data(interval(1),ii))/1e11,'.-r'); hold on;
        h2=stem(ind(ii),(FBCT_sel.data(interval(end),ii))/1e11,'.k'); hold on;
end

legend([h1,h2],datestr(xlim1,'HH:MM:SS'),datestr(xlim2,'HH:MM:SS'),'location','northwest'); legend boxoff


title(['Fill:',fill,' FBCT B',beam,' M=',num2str(Mbunches)])
xlim([1,3564])
xlabel('# bucket'); ylabel('N_b [10^{11} ppb] [um]');

if flagsave
    name=['FBCT_vs_buckets_zoom_B',beam,'V'];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='8';
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end


figure
set(gcf,'visible',flagshow);
interval_FBCT=find(FBCT_sel.time>xlim1 & FBCT_sel.time<xlim2);
interval_BSRT=find(BSRT_sel.time>xlim1 & BSRT_sel.time<xlim2);
for ii=1:length(ind)
        h1=stem(ind(ii),(FBCT_sel.data(interval_FBCT(1),ii))/1e11/(BSRT_sel.eH(interval_BSRT(1),ii)),'.-r'); hold on;
        h2=stem(ind(ii),(FBCT_sel.data(interval_FBCT(end),ii))/1e11/(BSRT_sel.eH(interval_BSRT(end),ii)),'.k'); hold on;
end

legend([h1,h2],datestr(xlim1,'HH:MM:SS'),datestr(xlim2,'HH:MM:SS'),'location','northwest'); legend boxoff


title(['Fill:',fill,' FBCT/BSRT B',beam,' M=',num2str(Mbunches)])
xlim([1,3564])
xlabel('# bucket'); ylabel('N_b [10^{11} ppb] / [um]');

if flagsave
    name=['BRILL_vs_buckets_zoom_B',beam,'V'];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='8';
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end



%% NAFF & FFT analysis
start_turns=1;
end_turns =BBQ.turns(end);
shift_turns=1e3;
window_turns=1e4; nharm=1;
BBQ=slideFFT(BBQ,start_turns,end_turns,shift_turns,window_turns,'FFT',nharm);
BBQ=slideFFT(BBQ,start_turns,end_turns,shift_turns,window_turns,'NAFF',nharm);
disp('Done!')
%% Plot  

flagsave=0;

showlines=0;
Nline=0;

if strcmp(plane,'V')
    xlim1=0.295;
    xlim2=0.325;
elseif strcmp(plane,'H')
    xlim1=0.275;
    xlim2=0.29;
else
    xlim1=0.2;
    xlim2=0.4;
end

f=figure(11);
h=imagesc(BBQ.FFT_freq(:,1),BBQ.FFT_slide_time(end:-1:1),BBQ.FFT_amp'); 
caxis([1 1e10])

ylim([BBQ.FFT_slide_time(1) BBQ.FFT_slide_time(end) ])
ylabel('Time'); xlabel('Q');
h=colorbar;
set(h,'yscale','log')
title(['Fill:',fill,' BBQ ',typeBBQ,' B',beam,plane])
hold on;
N=10;
ytime=linspace(BBQ.FFT_slide_time(1),BBQ.FFT_slide_time(end),N);
set(gca,'ytick',ytime)
set(gca,'yticklabel',datestr(ytime(end:-1:1),'HH:MM:SS'))



if isfield(BBQ,'NAFF_freq')
    for jj=1:nharm
        plot(BBQ.NAFF_freq(jj,end:-1:1),BBQ.NAFF_slide_time,'.k')
    end
end

figure(12)
ind1=[floor(size(BBQ.FFT_amp,2)/Nline):floor(size(BBQ.FFT_amp,2)/Nline):size(BBQ.FFT_amp,2)];
col_vec=distinguishable_colors(length(ind1));
for ii=1:length(ind1)
    figure(12)
    h=semilogy(BBQ.FFT_freq(:,1),BBQ.FFT_amp(:,ind1(ii)),'-','color',col_vec((ii),:)); hold on;
    if isfield(BBQ,'NAFF_freq')
        for jj=1:nharm
            plot(BBQ.NAFF_freq(jj,ind1(ii)),BBQ.NAFF_amp(jj,ind1(ii))*window_turns,'ok','markerfacecolor',col_vec((ii),:))
        end
    end
    
    
    if showlines
        figure(11); 
        hold on;
        h=hline(BBQ.FFT_slide_turns(ind1(ii)),'-');
        set(h,'color',col_vec(ii,:))
    end        
end

figure(11); xlim([xlim1 xlim2]); hold off;
figure(12); xlim([xlim1 xlim2]); hold off;

ylabel('BBQ amplitude [arb. units]'); xlabel('Q');
title(['Fill:',fill,'BBQ ',typeBBQ,' B',beam,plane])

figure(12)
if flagsave
    name=['spectrum_BBQ',typeBBQ,'_B',beam,plane];
    s=hgexport('readstyle','PRSTAB-10pt'); 
    hgexport(gcf,'',s,'applystyle',true);
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');  
end
%%
figure(11)
if ~flagsave
    name=['waterfall_BBQ',typeBBQ,'_B',beam,plane,'_',comment];
    s=hgexport('readstyle','PRSTAB-10pt'); 
    hgexport(gcf,'',s,'applystyle',true);
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');  
end

%% NAFF & Coll
flagsave=1;


NAFF.time=BBQ.NAFF_slide_time;
NAFF.data=BBQ.NAFF_freq;
NAFF.turns=BBQ.NAFF_slide_turns;

figure(1)
datall=NAFF.data(:);
NAFF.data2=[];
NAFF.time2=[];
    for jj=1:size(NAFF.data,2)
        for ii=1:size(NAFF.data,1)
        NAFF.data2=[NAFF.data2,NAFF.data(ii,jj)];
        
        end
        NAFF.time2=[NAFF.time2,NAFF.time(jj)+ 1e-6*rand(1,3)];
    end
    
[NAFF.time,ind]=sort(NAFF.time2);
NAFF.data=NAFF.data2(ind);

h=plot(NAFF.time,NAFF.data,'.k');

brush on;
pause
brush off;
x_remaining=get(h,'XData');
ind_del=~ismember(NAFF.time,x_remaining );
[~,ind]=find(ind_del==1);
NAFF.data(ind)=[];
NAFF.time(ind)=[];

h=plot(NAFF.time,NAFF.data,'-r');


ind_nan=isnan(NAFF.data);
NAFF.data(ind_nan)=[];
NAFF.time(ind_nan)=[];

[NAFF,Coll]=interp_timber_data(NAFF,Coll,1);



val=300;
NAFF.smooth_data=fastsmooth(NAFF.data,val,1,1);
Coll.smooth_data=fastsmooth(Coll.data,val,1,1);


figure()
h1=plot(NAFF.time,NAFF.data,'-','linewidth',1,'color',[.5 .5 .5]);
hold on;
h2=plot(NAFF.time,NAFF.smooth_data,'-k','linewidth',2);
h3=addaxis(Coll.time,Coll.data,'-r','linewidth',2);
h4=addaxisplot(Coll.time,Coll.smooth_data,2,'--r','linewidth',2);
addaxislabel(2,'Full gap [mm]')
addaxislabel(1,['Q',plane])

gmin=1.05*min(Coll.data);
gmax=0.95*max(Coll.data);

ind=find(Coll.smooth_data<gmin);
turn_sample_min=Coll.time(ind);
sm_tune_vec_min=NAFF.smooth_data(ind);

ind=find(Coll.smooth_data>gmax);
turn_sample_max=Coll.time(ind);
sm_tune_vec_max=NAFF.smooth_data(ind);

h5=vline(Coll.time(val),'-g');set(h5,'linewidth',1)
h5=vline(Coll.time(end-val),'-g');set(h5,'linewidth',1)

ind=find(turn_sample_min>Coll.time(val) & turn_sample_min<Coll.time(end-val));
turn_sample_min_fit=turn_sample_min(ind);
sm_tune_vec_min_fit=sm_tune_vec_min(ind);

ind=find(turn_sample_max>Coll.time(val) & turn_sample_max<Coll.time(end-val));
turn_sample_max_fit=turn_sample_max(ind);
sm_tune_vec_max_fit=sm_tune_vec_max(ind);

h6=addaxisplot(turn_sample_min_fit ,sm_tune_vec_min_fit ,1,'.g','linewidth',1);
h6=addaxisplot(turn_sample_max_fit,sm_tune_vec_max_fit,1,'.g','linewidth',1);

[fmin,~]=fit(turn_sample_min_fit' ,sm_tune_vec_min_fit','poly1');
f1=plot(turn_sample_min_fit',fmin.p1*turn_sample_min_fit'+fmin.p2,'--g','linewidth',2);
cmin=confint(fmin);

[fmax,~]=fit(turn_sample_max_fit' ,sm_tune_vec_max_fit','poly1');
f2=plot(turn_sample_max_fit',fmax.p1*turn_sample_max_fit'+fmax.p2,'--g','linewidth',2);
cmax=confint(fmax);

fitted=((fmin.p1)*(turn_sample_max_fit')+fmin.p2);
tsmax=(sm_tune_vec_max_fit')-fitted;
fitted=((fmax.p1)*(turn_sample_min_fit')+fmax.p2);
tsmin=fitted-(sm_tune_vec_min_fit');

ts=mean([tsmax;tsmin]);
ts_err=std([tsmax;tsmin]);

disp(['ts=',num2str(ts),' +/- ',num2str(ts_err)]);
ts_str=['ts=(',sprintf('%.1f',ts*1e4),' \pm ',sprintf('%.1f',ts_err*1e4),') 10^{-4}'];
disp(ts_str);
title(['B',beam,plane,' TCSG7, ',ts_str])

box on;

leg=legend([h1,h2,h3,h4,h5,h6,f1],'Q data',['Q ',num2str(val),'-smooth'],'Coll. data',['Coll ',num2str(val),'-smooth'],'Range','Q selected','Fit');
set(leg,'location','northwest')

N=10;
xtime=linspace(Coll.time(1),Coll.time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'))
set(gca,'xtick',(xtime))
rotateticklabel(gca,90,8);
xlim([Coll.time(1),Coll.time(end)]);

set(gcf,'position',[850   -98   964   384]); 

if flagsave
    name=['TCSG7_tuneshift_B',beam,plane];
    s=hgexport('readstyle','PRSTAB-6pt'); s.FixedFontSize='10'; 
    
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end


%% Risetimes in in TD from BBQ raw

close all;
flagsave=1;
flagshow='on';


figure(111)
plot(BBQ.time,BBQ.data); hold on
N=20;
xtime=linspace(BBQ.time(1),BBQ.time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'))
set(gca,'xtick',(xtime))
rotateticklabel(gca,90,8);
axisxx(BBQ.turns,BBQ.data);

%%
BBQ.rise_turns=1:1e6; % window on BBQ
BBQ.rise_time=BBQ.time(BBQ.rise_turns);
BBQ.rise_data=BBQ.data(BBQ.rise_turns);
figure(2)
set(gcf,'visible',flagshow)
up=envelope(BBQ.rise_turns,BBQ.rise_data,1000,'top');

h1=plot(BBQ.rise_time,BBQ.rise_data,'-k'); hold on;
h3=plot(BBQ.time,BBQ.data,'-'); set(h3,'color',[.5 .5 .5]);
h2=plot(BBQ.rise_time,up,'-r','linewidth',2); hold off;

time_2fit=BBQ.time(BBQ.rise_turns);
turns_2fit=1:length(BBQ.rise_turns);
L_all_2fit=up(turns_2fit);
ylabel('<x> [a.u.]')
title(['Fill:',fill,' B',beam,plane])
legend([h1,h2],'data','envelope')
N=10;
ytime=linspace(BBQ.FFT_slide_time(1),BBQ.FFT_slide_time(end),N);
set(gca,'xtick',ytime)
set(gca,'xticklabel',datestr(ytime,'HH:MM:SS'))
rotateticklabel(gca,90,6);

if flagsave
    name=['Envelope_instability'];
    s=hgexport('readstyle','PRSTAB-6pt');% s.FixedFontSize='20';s.width='20'; set(gcf,'PaperType','A3')
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end

%% fit with exponential
f=ezfit(turns_2fit'/machine.f0,L_all_2fit','y(x)=b+a*exp((x)/tau); b=1e6; tau=10; a=100');

figure(3)
set(gcf,'visible',flagshow)
h1=plot(turns_2fit,L_all_2fit,'--b','linewidth',3); hold on;
h2=plot(turns_2fit,f.m(2)+f.m(1)*exp((turns_2fit)/(f.m(3)*machine.f0)),'r','linewidth',3); 
h3=plot((turns_2fit),BBQ.rise_data(turns_2fit),':k','linewidth',1); hold off;

legend([h3,h1,h2],'Data','Envelope',['Exp. fit.: tau=',num2str(f.m(3))],'location','northwest');
legend boxoff;
xlabel('Turns');
ylabel('<x> [arb. units]')
title(['Fill:',fill,' B',beam,plane])

if flagsave
    name=['Risetime_',comment,'_B',beam,plane];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='15';
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end


%% Risetimes in in TD from  NAFF amplitudes

close all;
flagsave=1;
flagshow='on';

figure(111)
plot(BBQ.FFT_slide_time,BBQ.NAFF_amp); hold on
N=20;
xtime=linspace(BBQ.FFT_slide_time(1),BBQ.FFT_slide_time(end),N);
set(gca,'xticklabel',datestr(xtime,'HH:MM'))
set(gca,'xtick',(xtime))
rotateticklabel(gca,90,8);

%%

BBQ.rise_turns=find([BBQ.FFT_slide_turns>1 & BBQ.FFT_slide_turns<1.1e6]==1); % window on BBQ
BBQ.rise_time=BBQ.FFT_slide_time(BBQ.rise_turns);
BBQ.rise_data=BBQ.NAFF_amp(BBQ.rise_turns);
BBQ.rise_turns=BBQ.FFT_slide_turns(BBQ.rise_turns);

indnan=isnan(BBQ.rise_data);
BBQ.rise_data(indnan)=[];
BBQ.rise_turns(indnan)=[];
BBQ.rise_time(indnan)=[];

figure(2)
set(gcf,'visible',flagshow)
h1=plot(BBQ.rise_time,BBQ.rise_data,'-k'); hold on;

turns_2fit=(BBQ.rise_turns);
L_all_2fit=BBQ.rise_data;
ylabel('<x> [a.u.]')
title(['Fill:',fill,' B',beam,plane])
legend([h1],'data')
N=100;
ytime=linspace(BBQ.FFT_slide_time(1),BBQ.FFT_slide_time(end),N);
set(gca,'xtick',ytime)
set(gca,'xticklabel',datestr(ytime,'HH:MM:SS'))
rotateticklabel(gca,90,6);

if flagsave
    name=['Envelope_instability'];
    s=hgexport('readstyle','PRSTAB-6pt');% s.FixedFontSize='20';s.width='20'; set(gcf,'PaperType','A3')
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end

% fit with exponential
f=ezfit(turns_2fit'/machine.f0,L_all_2fit','y(x)=b+a*exp((x)/tau); b=1e6; tau=6; a=1e2');

figure(3)
set(gcf,'visible',flagshow)
h1=plot(turns_2fit,L_all_2fit,'--b','linewidth',3); hold on;
h2=plot(turns_2fit,f.m(2)+f.m(1)*exp((turns_2fit)/(f.m(3)*machine.f0)),'r','linewidth',3); 
h3=plot((turns_2fit),BBQ.rise_data,':k','linewidth',1); hold off;

legend([h3,h1,h2],'Data','Envelope',['Exp. fit.: tau=',num2str(f.m(3))],'location','northwest');
legend boxoff;
xlabel('Turns');
ylabel('<x> [arb. units]')
title(['Fill:',fill,' B',beam,plane])

if flagsave
    name=['Risetime_',comment,'_B',beam,plane];
    s=hgexport('readstyle','PRSTAB'); s.FixedFontSize='15';
    hgexport(gcf,'',s,'applystyle',true);
    saveas(gcf, [ResultDir,name,'.fig'],'fig');
    hgexport(gcf, [ResultDir,name,'.pdf'],s,'Format','pdf');
    hgexport(gcf, [ResultDir,name,'.png'],s,'Format','png');
end


%%

e1=2;
N1=0.86;

em=2;
Nm=1.1;

I1=150;

Im=I1*(em/e1)*(Nm/N1)
