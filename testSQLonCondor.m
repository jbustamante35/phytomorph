function [res] = testSQLonCondor()
    res = 1;
    mksqlite('open','./dataBase.db');
    mksqlite('CREATE TABLE flinks (id INTEGER PRIMARY KEY AUTOINCREMENT,key,time,data,type,target_key,source_key)');
    mksqlite('CREATE TABLE tlinks (id INTEGER PRIMARY KEY AUTOINCREMENT,key,time,linkType,parent)');
end

%{

    func = cFlow('testSQLonCondor');
    func.setMCRversion('v930');
    res = func('a')
    auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
    auth = auth{1};
    func.submitDag(auth,50,50);
%}