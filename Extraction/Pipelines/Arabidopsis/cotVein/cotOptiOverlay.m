function [stop] = cotOptiOverlay(optimValues,state,Z,dispImage,trans0,fixed,moving,dispi)
    stop = false;
    if ~isempty(Z)
        here = 1;
    end
    if dispi
        if isstruct(optimValues)
            trans = optimValues.bestx;
        else
            trans = optimValues;
        end


        imshow(dispImage,[]);
        szf = size(fixed)/2;
        szm = size(moving)/2;
        hold on;

        plot(trans0(5)+szf(2),trans0(4)+szf(1),'r*');
        plot(trans0(5)+trans(5)+szf(2),trans0(4)+trans(4)+szf(1),'g*');


        init_dB = bwboundaries(moving > .8);
        init_dB = init_dB{1};
        init_dB = bsxfun(@minus,init_dB,szm);
        init_dB = [init_dB ones(size(init_dB,1),1)];

        if isstruct(optimValues)
            swarm = optimValues.swarm;
            pidx = randperm(size(swarm,1));
            swarm = swarm(pidx,:);
            for r = 1:10
                T = buildTrans(swarm(r,:)+trans0);
                dB = (T*init_dB')';
                plot(dB(:,2)+szf(2),dB(:,1)+szf(1),'b')
            end
        end

        T = buildTrans(trans+trans0);
        dB = (T*init_dB')';
        plot(dB(:,2)+szf(2),dB(:,1)+szf(1),'g','LineWidth',3)

        %{
        for b = 1:size(curveBank,1)
            tmpTrans0 = curveBank(b,1:5);
            tmpTrans = curveBank(b,6:end);
            T = buildTrans(tmpTrans+tmpTrans0);
            dB = (T*init_dB')';
            plot(dB(:,2)+szf(2),dB(:,1)+szf(1),'r')
        end
        %}

        hold off
        drawnow





    end
end