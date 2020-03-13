function [] = cmdLineLogin()
    v = 1;
    while v ~= 0
        userName = input(['Please enter your user name.\n'],'s');
        fprintf(['Please enter your password.\n']);
        cmd = 'read PW;echo $PW';
        [~,pw] = system(cmd);
        pw = strtrim(pw);
        cmd = ['./iAuth_short ' userName ' ' pw];
        v = system(cmd);
        if v ~= 0
            fprintf(['Login not valid. Please try again.\n'])
        end
    end
end