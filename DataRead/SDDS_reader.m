function [data]=SDDS_reader(file_path)

data=[];
line_cou=0;
fillcou=1;
fid=fopen(file_path,'r');
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if regexp(tline,'&parameter') 
               newentry=regexp(tline,'name=(\w*),','tokens');
               name_newentry=newentry{1};
               item={};
               item.name=name_newentry;
               data=[data,item];
               line_cou=line_cou+1;
            end
            if regexp(tline,'&array') 
               newentry=regexp(tline,'name=(\w*),','tokens');
               name_newentry=newentry{1};
               item={};
               item.name=[name_newentry,'_length'];
               data=[data,item];
               item={};
               item.name=name_newentry;
               data=[data,item];
               line_cou=line_cou+1;
            end
            if regexp(tline,'&data') end;
            if regexp(tline,'!page') end;
            
            if isempty(regexp(tline,'&','match')) && isempty(regexp(tline,'SDDS','match')) &&  isempty(regexp(tline,'!','match'))
                if length(tline)> 150
                data(fillcou).data=strread(tline,'%f','delimiter',' ');
                else
                data(fillcou).data=tline;
                end
                fillcou=fillcou+1;
            end
           
        end
        fclose(fid);
end
