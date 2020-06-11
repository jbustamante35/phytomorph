function [m] = bwareaopenSmallHoles(bw,areaSmall,conn)
    if nargin == 2;conn = 4;end
    ibw = ~bw;
    fbw = bwareaopenSmall(ibw,areaSmall,conn);
    m = ~fbw;
end