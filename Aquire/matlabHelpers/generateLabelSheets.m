function [pdfFiles] = generateLabelSheets(labelCollection,messageTarget,csvFile,labelFormat,constFileData,constFolderData)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create object store
    global store;
    store = objectFSstore();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % list of output files
    pdfFiles = {};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % number of labels per sheet
    labelsPerSheet = labelFormat.numRows*labelFormat.numCols;
    % read in the text file for makinglabels
    labelData = readtext(csvFile);
    type = labelData(1,:);
    keys = labelData(2,:);
    labelData(1:2,:) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % remove later
    %textSize = round(mean(labelFormat.textBox(:,3:4)));
    %curdpi = 150;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make folder for the collection
    collectionLocation = [getenv('PHYTO_FILES') 'labelCollections' filesep];
    collectionLocation = [collectionLocation labelCollection.name filesep];
    mmkdir(collectionLocation);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % for each row in the csv file
    sheetCNT = 1;
    for v = 1:size(labelData,1)
        tm = clock;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make new sheet if needed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if mod(v-1,labelsPerSheet) == 0
            sheet = ones(labelFormat.imageSize);
            labelCNT = 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract the keys from the row
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s = struct;
        text = [];
        % loop over each key field
        for k = 1:numel(keys)
            s.(keys{k}) = labelData(v,k);
            % if human readable
            if contains(lower(type{k}),'h')
                text = [text keys{k} ':' labelData{v,k} '\n'];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make objects
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm2 = clock;
        % generate sample label
        sampleLabel = store.generate('label',text);
        % generate the sample message body
        sampleBody = store.generate('sampleBody',constFileData,constFolderData);
        % generate message for start tile - addressed to messageTarget
        sampleMsg = store.generate('daqMessage',sampleLabel,messageTarget,sampleBody);
        % generate envelope and add message
        sampleEnvelope = store.generate('envelope',messageTarget,sampleMsg);
        % add the envelope to the tile
        store.invoke(sampleLabel,'attachEnvelope',sampleEnvelope);
        % add the label to the collection
        store.invoke(labelCollection,'addLabel',sampleLabel);
        objectCreateTime = etime(clock,tm2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make the index for thd qr image
        cvec = round(labelFormat.qrBox(labelCNT,1):(labelFormat.qrBox(labelCNT,1)+labelFormat.qrBox(labelCNT,3)-1));
        rvec = round(labelFormat.qrBox(labelCNT,2):(labelFormat.qrBox(labelCNT,2)+labelFormat.qrBox(labelCNT,4)-1));
        shortL = min([numel(cvec) numel(rvec)]);
        cvec = cvec(1:shortL);
        rvec = rvec(1:shortL);
        % make the qr image
        qrImage = sampleLabel.qrImage([shortL shortL]);
        sheet(rvec,cvec) = rgb2gray(qrImage);
        % make the index for thd qr image
        cvec = round(labelFormat.textBox(labelCNT,1):(labelFormat.textBox(labelCNT,1)+labelFormat.textBox(labelCNT,3)-1));
        rvec = round(labelFormat.textBox(labelCNT,2):(labelFormat.textBox(labelCNT,2)+labelFormat.textBox(labelCNT,4)-1));
        % make the text
        textImage = makeLabelText(text,labelFormat.textBox(labelCNT,3:4),12);
        sheet(rvec,cvec) = textImage;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        
        labelCNT = labelCNT + 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save sheet if needed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if mod(v,labelsPerSheet) == 0
            sheetImageName = [collectionLocation 'sheet' num2str(sheetCNT) '.tif'];
            sheetPDFName = [collectionLocation 'sheet' num2str(sheetCNT) '.pdf'];
            imwrite(sheet,sheetImageName);
            CMD = ['convert ' sheetImageName ' ' sheetPDFName];
            system(CMD);
            pdfFiles{end+1} = sheetPDFName;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % timing and reporting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dtm = etime(clock,tm);
        fprintf(['Done with:' num2str(v) ':' num2str(size(labelData,1)) '\n']);
        fprintf([num2str(objectCreateTime) ':' num2str(dtm) '\n'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    
end

%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% processTemplate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FilePath = '/mnt/scratch1/phytomorph_dev/Aquire/pdfTemplates/';
    FileList = {};
    FileExt = {'pdf'};
    FileList = gdig(FilePath,FileList,FileExt,1);

    height = .9;
    width = @(height)height;
    dx = 30;
    qr_text_spacing = 10;
    clear templateData;
    for e = 1:numel(FileList)
        templateData(e) = processAveryTemplate(FileList{e},height,width,dx,qr_text_spacing);
    end

    generateLabelSheets('/mnt/scratch1/phytomorph_dev/Aquire/phytoMeta/test.csv',templateData(1),'','');

%}