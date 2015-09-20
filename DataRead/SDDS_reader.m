function [data]=SDDS_reader(file_path)

data=[];
line_cou=0;
fillcou=1;
disp(file_path)
fid=fopen(file_path,'r');
        while 1
            tline = fgetl(fid);
            
            if ~ischar(tline), break, end
            if regexp(tline,'&parameter') 
%                disp(tline)
               newentry=regexp(tline,'name=(\w*)','tokens');
               name_newentry=newentry{1};
               item={};
               item.name=name_newentry;
               data=[data,item];
               line_cou=line_cou+1;
            
            elseif regexp(tline,'&array') 
%                disp(tline);
               newentry=regexp(tline,'name=(\w*)','tokens');
               name_newentry=newentry{1};
               item={};
               item.name=[name_newentry];
               data=[data,item];
               line_cou=line_cou+1;
            
            elseif regexp(tline,'&data')
            elseif regexp(tline,'!page')
            
            elseif isempty(regexp(tline,'&','match')) && isempty(regexp(tline,'SDDS','match')) &&  isempty(regexp(tline,'!','match'))
                try
                    data(fillcou).data=strread(tline,'%f','delimiter',' ');
                catch
                    data(fillcou).data=tline;
                end
                fillcou=fillcou+1;
            end
           
        end
        fclose(fid);
end