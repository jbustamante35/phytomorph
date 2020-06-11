classdef bAffine < basisT
    
    properties
        stretch;
    end
    
    methods
        function [obj] = bAffine(T,stretch)
            obj@basisT(T);
            obj.stretch = stretch;
        end
        
        
        function [] = plot(this,mag)
            if nargin == 1;mag = 10;end
            mid = this.E(3,:);
            t = mag*this.E(1,:);
            n = mag*this.E(2,:);
            plot(mid(1),mid(2),'k.');
            quiver(mid(1),mid(2),t(1),t(2),'g');
            quiver(mid(1),mid(2),n(1),n(2),'r');
        end
    end
    
end