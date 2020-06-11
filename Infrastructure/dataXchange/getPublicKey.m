function [keyFile] = getPublicKey(name)
    mksqlite('open','/mnt/scratch1/phytomorph_dev/Infrastructure/dataXchange/Xchange.db');
    collaborator = mksqlite(['select * from collaborator where name=''' name ''';']);
    mksqlite('close');
    keyPath = ['/chtc/phytomorphservice/collaborator/' collaborator.userID '/keys/public.key'];
    keyFile = dig(keyPath,{},'key');
end