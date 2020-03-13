function [d] = myImageContrast(moving,fixed,d,trans,trans0,disp)
    if nargin <= 4
        trans0 = 0;
        disp = false;
    end

    if disp
        imshow(fixed,[]);
        szf = size(fixed)/2;
        szm = size(moving)/2;
        hold on;
        dB = bwboundaries(moving > .8);
        T = buildTrans(trans+trans0);
        dB = dB{1};
        dB = [dB ones(size(dB,1),1)];
        dB = (T*dB')';
        plot(trans(5)+trans0(5)+szf(2),trans(5)+trans0(5)+szf(2),'r*');
        plot(dB(:,2)+szf(2)-szm(2),dB(:,1)+szf(1)-szm(1),'r')
        hold off
        drawnow
    end


    itemplate = transformImage(moving,fixed,d,trans,trans0);
    
    %d = bsxfun(@minus,itemplate,moving);
    %d = squeeze(sum(sum(d.*d,1),2).^.5);
    d = norm(moving(:) - itemplate(:));
    %d = -moving(:)'*itemplate(:);
end