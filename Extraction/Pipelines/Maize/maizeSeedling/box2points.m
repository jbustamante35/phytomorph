function [points] = box2points(box)
    points = [box(1:2);box(1:2) + [box(3) 0];box(1:2) + box(3:4);box(1:2) + [0 box(4)];box(1:2)];
end