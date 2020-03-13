classdef usbPort < doid & matlab.mixin.SetGet
    properties
        portN;
        computer;
    end
    
    methods
        function [obj] = usbPort(portN)
            obj = obj@doid();
            if nargin == 1
                obj.portN = portN;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % attach this to a computer
        % default to the computer attachment interface
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = attachComputer(obj,computer)
            % use the interface from computer
            computer.attachUSB(obj);
        end
        
         
    end
end