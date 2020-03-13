function [J] = myJacobian(func,x,h)
    %h = h'.*(x.*eps(x).^.5)';
    %h = x.*h;
    %ux = max(abs([x;tx]),[],1); % use x
    %h = h.*ux.*eps(ux).^.5;
    %h = .00001*ux;
    parfor e = 1:numel(x)
        delta = zeros(size(x));
        delta(e) = h(e);
        J(e,:) = (func(x+delta) - func(x-delta))/2*h(e);
    end
end