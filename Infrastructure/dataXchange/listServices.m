function [servicelist] = listServices()
    mksqlite('open','/mnt/scratch1/phytomorph_dev/Infrastructure/dataXchange/Xchange.db');
    servicelist = mksqlite('select * from service');
    mksqlite('close');
end