function [ret] = checkName(name,list)
    f = contains(list,name);
    qrName = readtext(list{find(f)});
   
    ret = strcmp(qrName,name);
end