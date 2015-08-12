function [D1_new,D2_new]=interp_timber_data(D1,D2,interp_index)
% function [D1_new,D2_new]=interp_timber_data(D1,D2,interp_index)
% if interp_index=1 interpolates data in D2 over timing in D1 and
% viceversa.
    D1_new={}; D2_new={};
    
    
    if interp_index==1
        for ii=1:size(D2.data,2)
            D2_new.data(:,ii)=interp1(D2.turns,D2.data(:,ii),D1.turns,'linear','extrap');
            D2_new.turns=D1.turns;
            D1_new=D1;
        end
    elseif interp_index==2
        for ii=1:size(D1.data,2)
            D1_new.data(:,ii)=interp1(D1.turns,D1.data(:,ii),D2.turns,'linear','extrap');
            D1_new.turns=D2.turns;
            D2_new=D2;
        end
    end

end