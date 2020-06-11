function [P] = randomPDF(n,sidx)
    distributionList = makedist();
    
    %{
    distributionList(1).name = 'Beta';
    distributionList(1).ni = 2;
    
    distributionList(2).name = 'Binomial';
    distributionList(2).ni = 2;
    
    distributionList(3).name = 'BirnbaumSaunders';
    distributionList(3).ni = 2;
    
    distributionList(4).name = 'Normal';
    distributionList(4).ni = 2;
    distributionList(4).fixed{1} = 0;
    %}
    
    if nargin == 1; sidx = randi(numel(distributionList),1,n);end
    
   
    P = distributionList(sidx);
    d = makedist(P{1});
    
    for e = 1:numel(P)
       % P(e).prob = @(x,y)pdf(P(e).name,
    end
    
end
