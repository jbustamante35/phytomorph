function [GT] = sampleAlongContours(inFilePath)
    try
        %%% scan for new images
        FileList = {};
        FileExt = {'tiff','TIF','tif'};
        verbose = 1;
        SET = sdig(inFilePath,FileList,FileExt,verbose);


        %%% for sampling curves
        T = 30;
        R = 50;
        nhSZ = [T R];
        [THETA RAD] = ndgrid(linspace(-pi,pi,T),linspace(0,30,R));
        NH = [RAD(:).*cos(THETA(:)) RAD(:).*sin(THETA(:))]';

        segmentSize = 61;    
        particleStore = labeledCurveStore();
        moleculeStore = curveStore();

        GV = [];
        cntGT = 1;

        %%% analysis of single image
        for e = 1:numel(SET)


            image = maizeImage(SET{e}{1});
            kernelCenters = image.extractKernelCenterPoints();
            curvesSet = image.extractClosedContours(@(x)maizeContour(x),10,1,[100 500]);
            curvesSet = curveStore.containsPoints(curvesSet,kernelCenters,@(x)any(x));
            containmentMap = curveStore.containmentMap(curvesSet);
            repCurves = curveStore.generateRepCurve(containmentMap,curvesSet,5);

            curves{e} = repCurves;
            for c = 1:numel(curves{e})
                try
                    atomicCurves = curves{e}(c).atomize(NH,nhSZ,segmentSize);
                    labels = particleStore.predictTipParticle(atomicCurves);
                    ujLabels = particleStore.predictJunctionUpper(atomicCurves);
                    ljLabels = particleStore.predictJunctionLower(atomicCurves);
                    upperJ = find(ujLabels);
                    lowerJ = find(ljLabels);
                    fidx = find(labels);
                    if ~isempty(fidx)
                        curves{e}(c).tagTip(fidx);
                        curves{e}(c).tagJunctions([upperJ lowerJ]);
                        tmpI = curves{e}(c).getImage();
                        imshow(tmpI,[]);
                        offSet = curves{e}(c).getOffSet();
                        hold on
                        % particleStore.selfDistance()
                        plot(curves{e}(c).segment(1,upperJ)-offSet(1),curves{e}(c).segment(2,upperJ)-offSet(2),'b*');
                        plot(curves{e}(c).segment(1,lowerJ)-offSet(1),curves{e}(c).segment(2,lowerJ)-offSet(2),'b*');
                        plot(curves{e}(c).junction(1,:)-offSet(1),curves{e}(c).junction(2,:)-offSet(2),'g');
                        plot(curves{e}(c).segment(1,:)-offSet(1),curves{e}(c).segment(2,:)-offSet(2),'r');
                        plot(curves{e}(c).segment(1,fidx)-offSet(1),curves{e}(c).segment(2,fidx)-offSet(2),'g*');
                        plot(curves{e}(c).midline(1,:)-offSet(1),curves{e}(c).midline(2,:)-offSet(2),'r');
                        drawnow
                        waitforbuttonpress
                        hold off
                        
                        
                        
                        
                        
                    end
                    button = questdlg('Accept Labels?');
                    if strcmp(button,'No')
                        curves{e}(c).labelCurve();
                    end
                    
                    
                    atomicCurves = curves{e}(c).labelAtoms(atomicCurves);
                    particleStore.insertParticles(atomicCurves);
                    
                    
                    button = questdlg('Insert curve into store?');
                    if strcmp(button,'Yes')
                        moleculeStore.insertParticles(curves{e}(c));
                    end
                    
                    
                    
                    GT(cntGT) = curves{e}(c);
                    cntGT = cntGT + 1;
                    
                    close all
                    figure;
                    maizeContour.plotCurveSet(GT);
                    axis equal
                    waitforbuttonpress
                    
                    
                    %{
                   figure;
                    for i = 1:size(S,1)





                        raw = GV(i,:);
                        raw = reshape(raw,vecSZ);
                        tmp = reshape(S(i,:),vecSZ);
                        plot(tmp(1,:),tmp(2,:),'r');
                        axis equal
                        hold on
                        plot(raw(1,:),raw(2,:),'b');
                        waitforbuttonpress


                    end
                    %%
                    INIT = C(1,:);
                    FINAL = C(end,:);
                    for i = 1:size(INIT,2)
                        DX(:,i) = linspace(INIT(i),FINAL(i),10);
                    end

                    for i = 1:size(DX,1)
                        CURVE = PCA_BKPROJ(DX(i,:),E,U);
                        CURVE = reshape(CURVE,vecSZ);
                    end
                    %}
                    
                    B = maizeContour.decomposeCurveSet(GT);
                    f = curveStore.curveProb(GT(1:end-1),GT(end));
                    close all
                   
                    figure;plot(f);
                    %waitforbuttonpress;
                    
                catch ME
                        ME;
                end
            end






        end
    catch ME
        ME;
    end
end

function [curve] = extractCurves(fileName)
    Mask = createKernelEdgeMask(fileName);
    Mask = imfill(Mask,'holes');
    B = bwboundaries(Mask);
    I = double(imread(fileName));
    for e = 1:numel(B)
        curve(e) = maizeContour(flipud(B{e}'));
        curve(e).filename = fileName; 
    end
end
%{
    sampleAlongContours('/mnt/spaldingdata/Takeshi/allMaizeMovies/');
    sampleAlongContours('/mnt/spaldingdata/Takeshi/allMaizeMovies/agt mutant/');
    sampleAlongContours('/mnt/doane/doane/RIL_Data/999-sm-6-11-7-1000-lg-136-10-19/');
%}