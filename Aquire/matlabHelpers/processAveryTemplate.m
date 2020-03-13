function [templateData] = processAveryTemplate(pdfFile,height,width,dx,qr_text_spacing)
    % convert the pdf template to tiff and read
    [~,I] = pdf2tif(pdfFile);
    % threshold the image
    M = I == 255;
    % clear the border
    M = imclearborder(M);
    % get region props
    R = regionprops(M,'Image','BoundingBox','Centroid');
    % order the regions
    [R,numRow,numCol] = orderRegionProps(M,R);
    
    % make the label boxes and get the bounding box for each label
    per = .99;
    H = [];
    W = [];
    for e = 1:numel(R)
        X = sum(R(e).Image,1);
        Y = sum(R(e).Image,2);
        H(e,:) = [min(X(:)) max(X(:))];
        W(e,:) = [min(Y(:)) max(Y(:))];
        labelBox(e,:) = R(e).BoundingBox;
    end
    % get the mean of the dims for the labels
    uH = mean(H,1);
    uW = mean(W,1);
    % make the height and width based on means
    height = height*uH(1);
    width = width(height);

    for e = 1:numel(R)
        X = sum(R(e).Image,1);
        Y = sum(R(e).Image,2);
        fidx = X >= per*uH(1);
        fidy = Y >= per*uW(1);
        Yv = 1:numel(Y);
        Xv = 1:numel(X);
        cy = round(mean(Yv(fidy)));
        useableY = Yv((cy - height/2):(cy + height/2));
        useableX = Xv(fidx);

        maxWidth = useableX(end) - useableX(1);
       

        zX = zeros(size(X));
        zY = zeros(size(Y));

        fidx = useableX(1:width)+dx;
        fidy = useableY;

        qrBox(e,:) = [min(fidx) min(fidy) max(fidx)-min(fidx) max(fidy)-min(fidy)];
        qrBox(e,1:2) = qrBox(e,1:2) + labelBox(e,1:2);

        textWidth = maxWidth - qrBox(e,3);
        textWidth = textWidth - qr_text_spacing - dx*2;
        
        
        
        upperRight = qrBox(e,1:2) + [qrBox(e,3) 0] + [qr_text_spacing 0];

        textBox(e,:) = [upperRight textWidth qrBox(e,4)];
        
        %zX(fidx) = 1;
        %zY(fidy) = 1;
        %qrLabel = zY*zX;
        %imshow(qrLabel,[]);
        %drawnow
    end
    
   
    centroid = [];
    for e = 1:numel(R)
        centroid = [centroid;R(e).Centroid];
    end
    
    templateData.centroid = centroid;
    templateData.textBox = textBox;
    templateData.qrBox = qrBox;
    templateData.labelBox = labelBox;
    templateData.imageSize = size(I);
    templateData.height = height;
    templateData.width = width;
    templateData.numRows = numRow;
    templateData.numCols = numCol;
    viewAveryLabelTemplate(I,templateData);
    
end

function [nR,numRow,numCol] = orderRegionProps(M,R)
    M1 = sum(M,1);
    M1 = bindVec(M1);
    M1 = M1 > .9;
    mR1 = regionprops(M1);
    numCol = numel(mR1);
    numRow = numel(R)/numCol;
    cen = [];
    for e = 1:numel(R)
        cen = [cen;R(e).Centroid];
    end
    [~,sidx] = sort(cen(:,2));
    for r = 1:numRow
        g = sidx(1:numCol);
        tmpCen = [];
        for e = 1:numel(g)
            tmpCen = [tmpCen;R(g(e)).Centroid];
        end
        [~,gidx] = sort(tmpCen(:,1));
        g = g(gidx);
        
        str = (r-1)*numCol + 1;
        stp = str + numCol-1;
        
        
        nR(str:stp) = R(g);
        sidx(1:numCol) = [];
    end
end