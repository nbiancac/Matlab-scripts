function str=str_extract(string,before,after)
    str=regexp(string,before,'split');
    str=regexp(str(2),after,'split');
    str=char(str{1}(1));
end