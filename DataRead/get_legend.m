function legenda=get_legend(scenario,subscenario);

legenda=[];

for kk=1:length(subscenario)
    newstr=char(regexprep(subscenario(kk),repstr,''));
    newstr=regexprep([char(scenario(kk)),' ',newstr],'_',' ');
    lista.legend=[lista.legend,{newstr}];
end

end