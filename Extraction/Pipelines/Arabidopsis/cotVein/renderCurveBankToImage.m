function [] = renderCurveBankToImage(I,curveBank,template,toPath)

    imshow(I,[]);
    hold on
    szm = size(template)/2;
    szf = size(I)/2;
    % draw the curve bank
    init_dB = bwboundaries(template > .8);
    init_dB = init_dB{1};
    init_dB = bsxfun(@minus,init_dB,szm);
    init_dB = [init_dB ones(size(init_dB,1),1)];

    for bou = 1:size(curveBank,1)
        tmpTrans0 = curveBank(bou,1:5);
        tmpTrans = curveBank(bou,6:end);
        T = buildTrans(tmpTrans+tmpTrans0);
        dB = (T*init_dB')';
        plot(dB(:,2)+szf(2),dB(:,1)+szf(1),'r')

        text(mean(dB(:,2))+szf(2),mean(dB(:,1))+szf(1),num2str(bou),'BackgroundColor','w');
    end
    drawnow
    saveas(gca,[toPath 'cotMap.jpg']);
end