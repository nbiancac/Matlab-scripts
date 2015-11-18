function [titlestr,legenda]=get_title_legend(subscenario)

substr=regexp(subscenario,'_','split');
repstr=[];
for kk=1:length(subscenario)-1
    repstr=unique([intersect(substr{kk},substr{kk+1},'stable')],'stable');
end
titlestr=strjoin(repstr,' ');

legenda=[];

for kk=1:length(subscenario)
    newstr=char(regexprep(subscenario(kk),repstr,''));
    newstr=regexprep([char(scenario(kk)),' ',newstr],'_',' ');
    lista.legend=[lista.legend,{newstr}];
end

end