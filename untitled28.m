parpool(3)

%%

X = zeros(7,13);

   
spmd

 codist = codistributor1d(1,[],size(X));
    XX = codistributed(X,codist);
  
    %testCo([],labindex)
    [firstCol,lastCol] = globalIndices(codist,1);

    XX(firstCol,1) = 1;


    L = getLocalPart(XX);
    L(1,1) = 2;

    class(L)
iscodistributed(L)

iscodistributed(XX)
    %L(1,4) = 1;

%labindex
    D = codistributed.build(L,codist);
end

X = gather(XX);
%%
A = magic(4);   %replicated on all workers
D = codistributed(A, codistributor1d(1));
L = getLocalPart(D)
%%
 spmd
                    
 codist = codistributor1d(1,[],[inputStackHeight,1]);
                    labindex
    
                    inputArray = codistributed(dFluid.data,codist);
                    localData = getLocalPart(inputArray);
                    


                    localData

                    pFluid = digitalFluid('parallelFluid',localData);
                    pFluid.stopData.targetQueue = dFluid.stopData.targetQueue;
                    pFluid.toCompute = dFluid.toCompute;
                
                    inputStackHeight = size(pFluid);

end