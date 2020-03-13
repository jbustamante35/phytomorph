function [data_out] = coreSample(data_in,domain,affineSequence,domainSZ,disp)
    data_out = zeros([domainSZ size(affineSequence,3)]);
    for e = 1:size(affineSequence,3)
        tmpT = affineSequence(:,:,e);
        tmpDomain = mtimesx(tmpT,domain);
        
        
        %{
        tmIDX = tmpDomain(3,1);
        frame = data_in(:,:,tmIDX);
        tmpI = ba_interp2(frame,tmpDomain(1,:),tmpDomain(2,:));
        %}
       
        
        tmpI = ba_interp3(data_in,tmpDomain(1,:),tmpDomain(2,:),tmpDomain(3,:));
        
        
        tmpI = reshape(tmpI,domainSZ);
        
        tmpI = tmpI';
        
        
        data_out(:,:,e) = tmpI;
        
        
       
        
        if disp
            v1 = tmpT(1:2,1);
            v2 = tmpT(1:2,2);
            p = tmpT(1:2,4);
           
            
            reSZ = 3;
            dispI = imresize(data_out(:,:,e),reSZ);
            tmIDX = tmpDomain(3,1);
            frame = data_in(:,:,tmIDX);
            frame(1:size(dispI,1),1:size(dispI,2)) = dispI;
            imshow(frame,[0 1]);
            hold on
            mag = 20;
            quiver(p(1),p(2),mag*v1(1),mag*v1(2),'Color','r')
            quiver(p(1),p(2),mag*v2(1),mag*v2(2),'Color','b')
            
            %plot(tmpDomain(1,:),tmpDomain(2,:),'.')
            drawnow
            hold off
        end
        
    end
end