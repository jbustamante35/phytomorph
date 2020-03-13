function [cTrace] = traceVec(F,grd,dms,cTrace,bvec,dl,alpha,beta)
    for l = 1:numel(dl)
        dM = zeros(size(grd,5),size(grd,6));
        dMsz = size(dM);
        IJ = cartprod(1:size(grd,5),1:size(grd,6));
        for k = 1:size(IJ,1)
            dM(k) = F{IJ(k,1),IJ(k,2)}(cTrace(l,:));
        end
        dM = reshape(dM,dMsz);
        %{
        for m = 1:size(grd,5)
            parfor p = 1:size(grd,6)
                %{
                dM(m,p) = interpn(dms{1},dms{2},dms{3},dms{4},grd(:,:,:,:,m,p),...
                    cTrace(l,1),cTrace(l,2),cTrace(l,3),cTrace(l,4),'makima');
                %}
                
            end
        end
        %}
        %vec = (inv(dM)'*bvec*dl(l)*alpha)';
        vec = (inv(dM')*bvec*dl(l)*alpha)';
        n = norm(vec);
        if n~= 0
            vec = vec / norm(vec);
        end
        cTrace(l+1,:) = cTrace(l,:) + beta*vec;
        l;
    end
end