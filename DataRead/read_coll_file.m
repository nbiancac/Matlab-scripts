function [names1,gaps1,betx1,bety1,nsig1]=read_coll_file(filepath,type)

% reads collimator files for LHC and HLLHC and gives back name, gap and
% beta functions:
% [names1,gaps1,betx1,bety1]=read_coll_file(filepath)

    fid=fopen([filepath]);
    tline = fgetl(fid);
    tline_dummy = fgetl(fid);
    fclose(fid);
    header_vec=regexp(tline,'\t','split');
    betx_ind=find(ismember(header_vec,'Betax[m]'));
    bety_ind=find(ismember(header_vec,'Betay[m]'));
    name_ind=find(ismember(header_vec,'Name'));
    material_ind=find(strcmp(header_vec,'Material'));
    gap_ind=find(strcmp(header_vec,'Halfgap[m]'));
    header_vec
    index=(strfind(header_vec,'nsig'));
    nsig_ind=find(~cellfun(@isempty,index));
    
    
    
    form=[];
    for ii=1:length(header_vec);
        if ii==name_ind || ii==material_ind
            form=[form,'%s'];
        else
            form=[form,'%f'];
        end
    end
    
    fid=fopen([filepath]);
    L=textscan(fid,form,'headerlines',1); 
    fclose(fid);
    
    if ~isempty(nsig_ind)
        [names1,gaps1,betx1,bety1,nsig1]=deal(L{name_ind},L{gap_ind},L{betx_ind},L{bety_ind},L{nsig_ind});
    else
        [names1,gaps1,betx1,bety1,nsig1]=deal(L{name_ind},L{gap_ind},L{betx_ind},L{bety_ind},nan.*ones(length(L{betx_ind}),1));
    end
    
    [names1,ia2,ic2] = unique(names1,'stable');
    [gaps1,betx1,bety1,nsig1]=deal(gaps1(ia2),betx1(ia2),bety1(ia2),nsig1(ia2));

    if ~strcmp(type,'all')
        index=(strfind(names1,char(type)));
        ind_coll=find(~cellfun(@isempty,index));
        [names1,gaps1,betx1,bety1,nsig1]=deal(names1(ind_coll),gaps1(ind_coll),betx1(ind_coll),bety1(ind_coll),nsig1(ind_coll));
    end

end