function [id] = addService(name)
    mksqlite('open','/mnt/scratch1/phytomorph_dev/Infrastructure/dataXchange/Xchange.db');
    [~,serviceID] = system('uuidgen');serviceID = strtrim(serviceID);
  
    mksqlite(['INSERT OR IGNORE INTO service (service,serviceID) '...
        'VALUES (''' name ''',''' serviceID ''');']);
    
    id = mksqlite(['SELECT id FROM service WHERE serviceID=''' serviceID ''';']);
    id = id.id;
    mksqlite('close');
end