function [w] = b2q(v)
    w = [1 0];
    if v
        w = [0 1];
    end
    %{
    w = [0 1];
    if v
        w = [1 0];
    end
    %}
end