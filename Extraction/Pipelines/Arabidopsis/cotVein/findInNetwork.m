function [q1] = findInNetwork(networkList,source)
     q1 = find(networkList(:,1) == source(1) & networkList(:,2) == source(2));
end