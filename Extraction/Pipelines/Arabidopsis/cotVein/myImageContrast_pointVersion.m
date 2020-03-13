function [d] = myImageContrast_pointVersion(moving,fixed,d,trans,trans0,disp,edgePointsM)
    if nargin <= 4
        trans0 = 0;
        disp = false;
    end
    %{
    if disp
        imshow(fixed,[]);
        szf = size(fixed)/2;
        szm = size(moving)/2;
        hold on;
        dB = bwboundaries(moving > .8);
        T = buildTrans(trans+trans0);
        dB = dB{1};
        dB = bsxfun(@minus,dB,szm);
        dB = [dB ones(size(dB,1),1)];
        dB = (T*dB')';
        plot(trans0(5)+szf(2),trans0(4)+szf(1),'r*');
        plot(trans0(5)+trans(5)+szf(2),trans0(4)+trans(4)+szf(1),'g*');
        plot(dB(:,2)+szf(2),dB(:,1)+szf(1),'g')
        hold off
        drawnow
    end
    %}

    itemplate = transformImage(moving,fixed,d,trans,trans0);




    v = ba_interp2(itemplate,edgePointsM.Points(:,2),edgePointsM.Points(:,1));
    extra = sum(v < .4);
    %itemplate = imresize(itemplate,.35);
    [templateP(:,1),templateP(:,2)] = find(edge(itemplate));
    [~,delta1] = nearestNeighbor(edgePointsM,templateP);

    %tic
    %delta1 = pdist2(edgePointsM.Points,templateP,'euclidean','Smallest',1);
    %toc
    %d = mean(delta1);
    delta2 = pdist2(templateP,edgePointsM.Points,'euclidean','Smallest',1);
    d = mean(delta1) + mean(delta2) + extra;
    %d = bsxfun(@minus,itemplate,moving);
    %d = squeeze(sum(sum(d.*d,1),2).^.5);
    %d = norm(moving(:) - itemplate(:));
    %d = -moving(:)'*itemplate(:);
end