q = parallel.pool.DataQueue;
W = 4;
pauseAmount = [5,1,1,1];
cmdValue = [1 0 1 0];

setBit = @(store,key,value)array(key) = value;

%afterEach(q,@(X)fprintf(X));

parfor e = 1:W
    cmdTester(q,cmdValue(e),pauseAmount(e))
end
%%
parpool(W);

spmd
    cmdValue(labindex)
end