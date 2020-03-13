function [] = recycleObjectArray(array)
    for e = 1:numel(array)
        recycleObject(array(e));
    end
end