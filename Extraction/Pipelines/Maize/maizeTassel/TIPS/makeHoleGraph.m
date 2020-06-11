function [holeGraph] = makeHoleGraph(skel)

    initSkel = skel;
    stop = false;
    N = 500;
    
    while ~stop
        nextSkel = initSkel;
        for e = 1:N
            endPoints = bwmorph(nextSkel,'endPoints');
            nextSkel = (endPoints == 0) & (nextSkel == 1);
        end
        
        nextSkel = bwmorph(nextSkel,'skeleton',inf);
        if all((nextSkel == initSkel))
            stop = true;
        end

        initSkel = nextSkel;
    end
    holeGraph = makeTasselGraph(initSkel);
end