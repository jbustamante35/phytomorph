function [] = xferp(sourceD,targetD)
    mmkdir(targetD);

    cmd = ['ils ' sourceD];
    [r,out] = system(cmd);
    out = strsplit(out,'\n');
    rm = [];
    for e = 1:numel(out)
        if ~isempty(out{e})
            if strcmp(out{e}(1),'/')
                rm = [rm e];
            else
                out{e} = [sourceD strtrim(out{e})];
            end
        else
            rm = [rm e];
        end
    end
    out(rm) = [];
    
    parfor e = 1:numel(out)
        tic
        cmd = ['iget ' out{e} ' ' targetD];
        res(e) = system(cmd);
        toc
    end
    
end

%{

    sourceD = '/iplant/home/guoshengwu/maizeData/coleoptileEmergence/mapping/20180605Camera6/';
    targetD = '/mnt/snapper/nate/redCapTest/20180605Camera6_ver2/';
    xferp(sourceD,targetD);


    

%}

