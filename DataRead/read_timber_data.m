function [time,data]=read_timber_data(variable,time_format)
% function [time,data]=read_timber_data(variable,time_format)
    disp(['Loading ',variable])
    if regexp(variable,'.csv')
        fid=fopen([variable]);
        L=textscan(fid,'%s\n','headerlines',3);
        fclose(fid);
        data=dlmread(variable,',',1,2);;
        
        time=[];
        for ii=1:length(L{1})
                time=[time,datenum(char(L{1}(ii)),time_format)];
          
        end
    elseif regexp(variable,'.dat')
        data=dlmread([variable],',',0,1);
        fid = fopen([variable]);
        C = textscan(fid,'%s%*[^\n]','delimiter',',');
        fclose(fid);
        time=[];
        for ii=1:length(C{1})
            time=[time,datenum(C{1}{ii},time_format)];
        end
    end
end