function [] = plotFrame(frame,mag,CL)
    if nargin == 2
        CL = {'r' 'b'};
    end
    quiver(frame(1,3),frame(2,3),frame(1,1),frame(2,1),mag,CL{1})
    quiver(frame(1,3),frame(2,3),frame(1,2),frame(2,2),mag,CL{2})
    plot(frame(1,3),frame(2,3),'g.')
end