function [pmd] = generatePhytomorphDictionary()

    pmd = struct;


    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % version and core delimiter values
    %%%%%%%%%%%%%%%%%%%%%%%%
    pmd.version = 1;
    
    pmd.objects.keys = {'uuid','genDate','type'};
    
    pmd.objects.type.values = {'msg','tile','triggerTile','file','folder'};
    
    pmd.objects.msg.keys = {'from','to','type','body'};
    
    

    %{
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % keys for headers
    %%%%%%%%%%%%%%%%%%%%%%%%
    pmd.header.keys.from = 'from';
    pmd.header.keys.to = 'to';
    pmd.header.keys.body = 'body';
    pmd.header.keys.type = 'type';
    
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%
    % keys for communicating where to save data on the file system
    %%%%%%%%%%%%%%%%%%%%%%%%
    pmd.fileSystem.keys.file = 'file';
    pmd.fileSystem.keys.folder = 'fldr';
   

    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % keys communicating about objects in the phytomorphDAQ
    % solution.
    %%%%%%%%%%%%%%%%%%%%%%%%
    % all objects have these available
    pmd.objects.keys.uuid = 'uuid';
    pmd.objects.keys.generateDate = 'genDate';
    pmd.objects.keys.type = 'otype';
    pmd.objects.values.type.tile = 'tile';
    pmd.objects.values.type.tileStrip = 'tileStrip';
    pmd.objects.values.type.camera = 'camera';
    pmd.objects.values.type.scanner = 'scanner';
    pmd.objects.values.type.computer = 'computer';
    
    
    
    % data for a tile
    pmd.objects.tileSequence.readme = ['keys and values for a qr tile'...
        'including sequenceID'];
    pmd.objects.tile.keys.sequenceID = 'seqID'; % link to parent sequence
    pmd.objects.tile.keys.type = 'tileType';
    
    
    pmd.objects.tile.values.type.qrstart = 'strMsg';
    pmd.objects.tile.values.type.qrend = 'stpMsg';
    pmd.objects.tile.values.type.metadata = 'metaMsg';
    pmd.objects.tile.values.type.sampledata = 'smplMsg';
    
    
    % data for a sequence of tiles
    pmd.objects.tileSequence.readme = ['keys and values for a sequence of qr tiles'...
        'including title, owner, date-generated, etc'];
    pmd.objects.tileSequence.keys.title = 'tsTitle';
    pmd.objects.tileSequence.keys.owner = 'tsOwner';
    pmd.objects.tileSequence.keys.numCapTiles = 'numCapTiles';
    pmd.objects.tileSequence.keys.numMetaTiles = 'numMetaTiles';
    pmd.objects.tileSequence.keys.numSampleTiles = 'numSampleTiles';
    
    
    
    pmd.header.keys.message ='msg';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % keys for device communication
    %%%%%%%%%%%%%%%%%%%%%%%%
    pmd.device.readme = ['key word dictionary for communcation from QR' ...
    'codes to devices on computer'];
    pmd.device.keys.usbPort = 'usbPort';
    pmd.device.keys.type = 'deviceType';

    %%%%%%%%%%%%%%%%%%%%%%%%
    % keys for tiggering the type of hash to use for file 
    % signing and verification
    %%%%%%%%%%%%%%%%%%%%%%%%
    pmd.device.readme = ['key word dictionary triggering the type of file' ...
        'verification, encryption, and string hashing used.'];
    pmd.verification.keys.defaultFile = 'sha512File';
    pmd.verification.keys.defaultString = 'sha512String';
    pmd.verification.keys.md5File = 'md5';
    pmd.verification.keys.sha256File = 'sha256';
    pmd.verification.keys.sha512File = 'sha512';
    pmd.verification.keys.md5String = 's-md5';
    pmd.verification.keys.sha256String = 's-sha256';
    
    %}
    
    %{
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % version and core delimiter values
    %%%%%%%%%%%%%%%%%%%%%%%%
    pmd.version = 1;
    pmd.delimiter.value.kv = '_';
    pmd.delimiter.value.json = ':';
    pmd.delimiter.value.sep = ',';
    pmd.delimiter.value.cellOpening = '{';
    pmd.delimiter.value.cellClosing = '}';
    pmd.delimiter.value.arrayOpening = '[';
    pmd.delimiter.value.arrayClosing = ']';

    %%%%%%%%%%%%%%%%%%%%%%%%
    % keys for headers
    %%%%%%%%%%%%%%%%%%%%%%%%
    pmd.header.keys.from = 'from';
    pmd.header.keys.to = 'to';
    pmd.header.keys.body = 'body';
    pmd.header.keys.type = 'type';
    
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%
    % keys for communicating where to save data on the file system
    %%%%%%%%%%%%%%%%%%%%%%%%
    pmd.fileSystem.keys.file = 'file';
    pmd.fileSystem.keys.folder = 'fldr';
   

    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % keys communicating about objects in the phytomorphDAQ
    % solution.
    %%%%%%%%%%%%%%%%%%%%%%%%
    % all objects have these available
    pmd.objects.keys.uuid = 'uuid';
    pmd.objects.keys.generateDate = 'genDate';
    pmd.objects.keys.type = 'otype';
    pmd.objects.values.type.tile = 'tile';
    pmd.objects.values.type.tileStrip = 'tileStrip';
    pmd.objects.values.type.camera = 'camera';
    pmd.objects.values.type.scanner = 'scanner';
    pmd.objects.values.type.computer = 'computer';
    
    
    
    % data for a tile
    pmd.objects.tileSequence.readme = ['keys and values for a qr tile'...
        'including sequenceID'];
    pmd.objects.tile.keys.sequenceID = 'seqID'; % link to parent sequence
    pmd.objects.tile.keys.type = 'tileType';
    
    
    pmd.objects.tile.values.type.qrstart = 'strMsg';
    pmd.objects.tile.values.type.qrend = 'stpMsg';
    pmd.objects.tile.values.type.metadata = 'metaMsg';
    pmd.objects.tile.values.type.sampledata = 'smplMsg';
    
    
    % data for a sequence of tiles
    pmd.objects.tileSequence.readme = ['keys and values for a sequence of qr tiles'...
        'including title, owner, date-generated, etc'];
    pmd.objects.tileSequence.keys.title = 'tsTitle';
    pmd.objects.tileSequence.keys.owner = 'tsOwner';
    pmd.objects.tileSequence.keys.numCapTiles = 'numCapTiles';
    pmd.objects.tileSequence.keys.numMetaTiles = 'numMetaTiles';
    pmd.objects.tileSequence.keys.numSampleTiles = 'numSampleTiles';
    
    
    
    pmd.header.keys.message ='msg';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % keys for device communication
    %%%%%%%%%%%%%%%%%%%%%%%%
    pmd.device.readme = ['key word dictionary for communcation from QR' ...
    'codes to devices on computer'];
    pmd.device.keys.usbPort = 'usbPort';
    pmd.device.keys.type = 'deviceType';

    %%%%%%%%%%%%%%%%%%%%%%%%
    % keys for tiggering the type of hash to use for file 
    % signing and verification
    %%%%%%%%%%%%%%%%%%%%%%%%
    pmd.device.readme = ['key word dictionary triggering the type of file' ...
        'verification, encryption, and string hashing used.'];
    pmd.verification.keys.defaultFile = 'sha512File';
    pmd.verification.keys.defaultString = 'sha512String';
    pmd.verification.keys.md5File = 'md5';
    pmd.verification.keys.sha256File = 'sha256';
    pmd.verification.keys.sha512File = 'sha512';
    pmd.verification.keys.md5String = 's-md5';
    pmd.verification.keys.sha256String = 's-sha256';
    pmd.verification.keys.sha512String = 's-sha512';
    %}
end

