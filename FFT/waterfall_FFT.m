disp(' ')
disp(' waterfall_FFT.m ')
disp(' ver 2.8  October 18, 2012 ')
disp(' by Tom Irvine  Email: tomirvine@aol.com ')
disp(' ')
disp(' This program calculates the one-sided, full amplitude FFT ')
disp(' of a time history ')
disp(' ')
disp(' The time history must be in a two-column matrix format: ')
disp(' Time(sec)  & amplitude ')
disp(' ')
%
close all;
%
clear THM;
clear tmx;
clear tmi;
clear ts;
clear te;
clear n;
clear n1;
clear n2;
clear dt;
clear po;
clear as;
clear ts;
clear amp;
clear tim;
clear aa;
clear bb;
clear Y;
clear store;
clear freq_p;
clear time_a;
clear store_p;
%
fig_num=1;
%
[t,y,dt,sr,tmx,tmi,n,ncontinue]=enter_time_history();
%
THM=[t y];
%
if(ncontinue==1)
%
    [tim,amp]=waterfall_FFT_time_process(tmi,tmx,dt,THM,n);
%
    clear THM;
%
    for ijk=1:10
%
        clear store;
        clear freq;
        clear store_p;
        clear freq_p;
        clear time_a;
        clear dt;
        clear df;
        clear NW;
        clear mmm;
%
        [dt,df,mmm,NW,io,minf,maxf]=waterfall_FFT_advise(tim,amp);
%
        [mk,freq,time_a,dt,NW]=...
                     waterfall_FFT_time_freq_set(mmm,NW,dt,df,maxf,tmi,io);
%
        [store,store_p,freq_p,max_a,max_f]=...
                            waterfall_FFT_core(NW,mmm,mk,freq,amp,minf,io);    
%
        [fig_num]=...
         waterfall_FFT_plots(NW,freq_p,time_a,store_p,max_f,max_a,fig_num);
%
        disp(' ')
        disp(' Repeat analysis with new parameters? ')
        repeat = input(' 1=yes  2=no  ');
%
        if(repeat==1)
           close all hidden;
           fig_num=1;
        else    
           break
        end
    end
%
    tim=fix_size(tim);
    amp=fix_size(amp);
%
    signal=[tim amp];
%
    out5 = sprintf('\n The time history is renamed as "signal" \n');
    disp(out5);      
%
end