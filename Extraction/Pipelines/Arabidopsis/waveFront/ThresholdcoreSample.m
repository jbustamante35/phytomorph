function [data_in] = ThresholdcoreSample(data_in)
    cp = ceil(size(data_in)/2);
    cp(end) = [];
    for e = 1:size(data_in,3)
        data_in(:,:,e) = data_in(:,:,e) > data_in(cp(1),cp(2),e);
    end
end