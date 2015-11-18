function [titlestr,legenda]=get_title_legend(scenario,subscenario)

substr=regexp(subscenario,'_','split');
%%
repstr=intersect(substr{1},substr{2},'stable');
for kk=2:length(subscenario)
    repstr=unique([intersect(substr{kk},repstr,'stable')],'stable');
end
titlestr=strjoin(repstr,' ');
%%
legenda=[];

for kk=1:length(subscenario)
    newstr=char(regexprep(subscenario(kk),repstr,''));
    newstr=regexprep([char(scenario(kk)),' ',newstr],'_',' ');
    legenda=[legenda,{newstr}];
end

end