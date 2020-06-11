function [id] = addXchangePoint(xchangeName,ownerName)
    [~,siteID] = system('uuidgen');siteID = strtrim(siteID);
    mksqlite('open','/mnt/scratch1/phytomorph_dev/Infrastructure/dataXchange/Xchange.db');
    ownerID = mksqlite(['SELECT id FROM collaborator WHERE name =''' ownerName ''';']);
    mksqlite(['INSERT OR IGNORE INTO xchangePoint (name,siteID,owner) VALUES ' ...
         '(''' xchangeName ''',''' siteID ''',''' num2str(ownerID.id) ''');']);
    id = mksqlite(['SELECT id FROM xchangePoint WHERE siteID=''' siteID ''';']);
    id = id.id;
    mksqlite('close');
end