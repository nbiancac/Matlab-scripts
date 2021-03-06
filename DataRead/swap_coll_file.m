function swap_coll_file(filepath,newfilepath,header2swap,values2swap);

% reads collimator files for LHC and HLLHC and gives back name, gap and
% beta functions:
% [names1,gaps1,betx1,bety1]=read_coll_file(filepath)

    fid=fopen([filepath]);
    tline = fgetl(fid);
    tline_dummy = fgetl(fid);
    fclose(fid);
    full_header=tline;
    header_vec=regexp(full_header,'\t','split');
    ind2swap=find(ismember(header_vec,header2swap));
    
    ind_vec=[];
    form=[];
    for ii=1:length(header_vec)
        ind_vec=[ind_vec,find(ismember(header_vec,char(header_vec(ii))))];
        
        if strcmp(header_vec(ii),'Name') || strcmp(header_vec(ii),'Material')
            form=[form,'%s'];
        else
            form=[form,'%f'];
        end
    end
    
    fid=fopen([filepath]);
    L=textscan(fid,form,'headerlines',1); 
    fclose(fid);
    
    L(ind2swap)={values2swap'};
    
    fid=fopen([newfilepath],'w');
    fprintf(fid,'%s\n',full_header);
    for ii=1:length(L{1})
        
        ev_str=[];
        for jj=1:length(header_vec)
            
            if strcmp(header_vec(jj),'Name') || strcmp(header_vec(jj),'Material')
                ev_str=[ev_str,[char(L{jj}(ii)),'   ']];
            else
                ev_str=[ev_str,[num2str(L{jj}(ii)),'   ']];
            end
            
        end
        
        fprintf(fid,'%s\n',(ev_str));
    end
    fclose(fid);
    
    disp('New file written: ')
    disp([newfilepath])
end