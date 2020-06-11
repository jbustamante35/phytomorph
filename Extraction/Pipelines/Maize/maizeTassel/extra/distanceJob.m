function [ret] = distanceJob(sourceData,targetData,domainData)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
domainW = domainData.xData(1);
domainW_num = domainData.xData(2);
domainL = domainData.yData(1);
domainL_num = domainData.yData(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make sample domain
[x1,x2] = ndgrid(linspace(-domainW,domainW,domainW_num),linspace(-domainL,domainL,domainL_num));
szX = size(x1);
x = [x1(:) x2(:) ones(size(x1(:)))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the search space function to map into affine transformations    
dX_range = [[50;0],[50;0]];
dR_range = [inf,0];
dS_range = [[2;.5],[2;.5]];
funcT = fwdT.makeTF(dX_range,dR_range,dS_range);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the function to create the frame at P
fT = @(P,gD)simpleAffine(P,gD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


target_rgbGlob = double(imread(targetData.rgb_fileName))/255;
target_maskGlob = double(imread(targetData.mask_fileName));


source_rgbGlob = double(imread(sourceData.rgb_fileName))/255;
source_maskGlob = double(imread(sourceData.mask_fileName));


for target = 1:size(targetData.pixelList,1)
   
   
    
    % set the target point
    targetPoint = targetData.pixelList(target,:);
    % make the target point object
    targetP = funcP(targetPoint,fT);
    % init the transformation
    targetP.getT(target_rgbGlob);
    % get the "function" - aka image patch at P
    targetF = targetP.getF(target_rgbGlob,x,szX);
    % get the mask at P
    targetM = targetP.getF(target_maskGlob,x,szX);
    
    
    for source = 1:size(sourceData.pixelList,1)
        
        % make a crop box for the source
        centerPoint = sourceData.pixelList(source,:);
        boxSZ =[domainW*2+100 domainL*2+100];
        upperLeft = centerPoint - boxSZ/2;
        
        % be careful - if we start to match rectangles
        bbox = [upperLeft boxSZ];
        
        box_rgbGlob = imcrop(source_rgbGlob,bbox);
        box_maskGlob = imcrop(source_maskGlob,bbox);
        
        
        % set the source point
        sourcePoint = hsize(box_maskGlob);
        % make the source point object
        sourceP = funcP(sourcePoint,fT);
        % init the transformation
        sourceP.getT(box_rgbGlob);
        
        
        % the zero point
        dT = [0 0 .1 1/3 1/3];
        % init the contrast function
        df = @(t)contrast1(t,funcT,sourceP,box_rgbGlob,box_maskGlob,targetF,targetM,x,szX);
        order = 4;
        ops = optimset('Display','iter','TolX',10^-order,'TolFun',10^-order);
        [tp(target,source,:),c(target,source)] = fminsearch(df,dT,ops);
        %{
        options = optimoptions('particleswarm','Display','iter','UseParallel',true,'SwarmSize',20);
        x = particleswarm(df,numel(dT),-1*ones(numel(dT),1),1*ones(numel(dT),1),options);
        %}
    end

    ret.tp = tp;
    ret.c = c;
end
