function [Z,midx] = spaceBLOCK(tmp,gz,Eq,Uq,cn2)
    % make the grid
    [n1,n2] = ndgrid(-gz:gz,-gz:gz);
    n = [n1(:) n2(:) ones(size(n1(:)))];

    % make mask
    scanSliceMask = (tmp(:,:,1) ~= 0);
    % rotate and remove the edges
    for e = 1:4
        scanSliceMask(:,1) = 0;
        scanSliceMask = imrotate(scanSliceMask,90);
    end
    % errode the mask
    scanSliceMask = imerode(scanSliceMask,strel('disk',11,0));
    % get the gradient
    [g1,g2] = gradient(tmp);
    % find the mask locations
    [m1,m2] = find(scanSliceMask);
    midx = find(scanSliceMask);
    MID = [m2 m1];
    % interpolate the grad
    G1 = ba_interp2(g1,m2,m1);
    G2 = ba_interp2(g2,m2,m1);
    % make the TAN and NOR
    TAN = [G1 G2];
    TAN = bsxfun(@times,TAN,sum(TAN.*TAN,2).^-.5);
    NOR = [TAN(:,2) -TAN(:,1)];
    % maket affine transformation
    affine = [];
    for e = 1:numel(m1)
        affine(:,:,e) = [[TAN(e,:)' NOR(e,:)' MID(e,:)'];[0 0 1]];
    end
    % permute affine and move domain(s)
    affine = permute(affine,[3 2 1]);
    nt = mtimesx(affine,n,'T');
    % remove the affine co-ordinate
    nt = nt(:,:,1:2);
    nt = permute(nt,[2 1 3]);
    szN = size(nt);

    nt = reshape(nt,[prod(szN(1:2)) szN(3)]);
    F = ba_interp2(tmp,nt(:,1),nt(:,2));
    F = reshape(F,[szN(1) szN(2)]);


    Fc = PCA_REPROJ_T(F,Eq,Uq,0);

    Z = zeros([size(tmp) cn2]);
    for k = 1:size(Z,3)
        tmpK = Z(:,:,k);
        tmpK(midx) = Fc(k,:);
        Z(:,:,k) = tmpK;
    end


end