function [time,data]=read_timber_data(variable,time_format)
% function [time,data]=read_timber_data(variable,time_format)
    data=dlmread([variable],',',0,1);
    fid = fopen([variable]);
    C = textscan(fid,'%s%*[^\n]','delimiter',',');
    fclose(fid);
    time=[];
    for ii=1:length(C{1})
    time=[time,datenum(C{1}{ii},time_format)];
    end

end