function [skeleton] = removeBranch(skeleton,networkPoints,rPath,branchPoints)
    % skeleton := image to remove points from
    % networkPoints := array of x-y ppoints
    % rPath := idx list of points to remove
    % branchPoints := [x-y] points to keep

    % get the path from the network
    PTH = networkPoints(rPath,:);
    % remove the branch points from the path
    PTH = setdiff(PTH,branchPoints,'rows');
    % remove the points from the skeleton image
    for p = 1:size(PTH,1)
        skeleton(PTH(p,1),PTH(p,2)) = 0;
    end
    skeleton = bwmorph(skeleton,'skeleton',inf);
end