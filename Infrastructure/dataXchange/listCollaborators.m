function [collaboratorlist] = listCollaborators()
    mksqlite('open','/mnt/scratch1/phytomorph_dev/Infrastructure/dataXchange/Xchange.db');
    collaboratorlist = mksqlite('select * from collaborator');
    mksqlite('close');
end