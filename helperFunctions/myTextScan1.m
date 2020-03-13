function [col] = myTextScan1(str)
    col = textscan(str,'%s','Delimiter','\t');
    col = col{1}';
end