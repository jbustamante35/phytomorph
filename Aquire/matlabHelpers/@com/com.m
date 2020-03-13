classdef com < doid
    
    properties
        hName;
        usbPorts;
    end
    
    methods
        function [obj] = com(hName)
            obj = obj@doid();
            if nargin == 1
                obj.hName = hName;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this function attaches USB port to computer
        % note that if a USB object is asked to attach
        % to computer it uses this interface
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = attachUSB(obj,usbObject)
            % transform pointer for usb
            [ptr] = dptr.transformInput(usbObject);
            % transform pointer for computer
            [computerObj] = dptr.transformInput(obj);
            % add the usbPort pointer to the array of usbports
            obj.usbPorts = [obj.usbPorts ptr];
            % add computer to port
            usbObject.computer = computerObj;
        end
        
        
        
    end
end