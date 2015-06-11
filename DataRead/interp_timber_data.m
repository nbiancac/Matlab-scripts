function [D1_new,D2_new]=interp_timber_data(D1,D2,interp_index)
% function [D1_new,D2_new]=interp_timber_data(D1,D2,interp_index)
% if interp_index=1 interpolates data in D2 over timing in D1 and
% viceversa.
    D1_new={}; D2_new={};
    if interp_index==1
        D2_new.data=interp1(D2.turns,D2.data,D1.turns);
        D2_new.turns=D1.turns;
    elseif interp_index==2
        D1_new.data=interp1(D1.turns,D1.data,D2.turns);
        D1_new.turns=D2.turns;
    end

end