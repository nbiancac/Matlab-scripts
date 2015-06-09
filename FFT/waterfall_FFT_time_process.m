%
%  waterfall_FFT_time_process.m  ver 1.0  October 18, 2012
%
function[tim,amp]=waterfall_FFT_time_process(tmi,tmx,dt,THM,n)
%
out3 = sprintf('\n  start  = %g sec    end = %g sec \n',tmi,tmx);
disp(out3)
disp(' Enter processing option for signal duration');
disp(' 1=whole time history   2=segment ');
po=input('');
%
if(po==2)
        disp(' ');
        ts=input(' Enter segment start time  ');
        te=input(' Enter segment end   time  ');
        disp(' ');
%
        if(ts>tmx)
            ts=tmx;
        end
        if(te<tmi)
            ts=tmi;
        end
%
        n1=fix((ts-tmi)/dt);
        n2=fix((te-tmi)/dt);
%
        if(n1<1)
            n1=1;
        end
%
        if(n2>n)
            n2=n;
        end
%
        if(n1>n2)
            n2=n1;
        end
%
        as=THM(:,2);
%
        ts=THM(:,1);
%
        amp=as(n1:n2)';
        tim=ts(n1:n2)';
%
else
        amp=THM(:,2);
        tim=THM(:,1);
end