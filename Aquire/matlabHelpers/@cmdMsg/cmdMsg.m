classdef cmdMsg < doid
    
    properties
        cmdType;
        cmd;
        varargin;
    end
    
    methods
        function [obj] = cmdMsg(cmdType,cmd,varargin)
            obj.cmdType = cmdType;
            obj.cmd = cmd;
            obj.varargin = varargin;
        end
        
        function [] = buildCmd(obj)
            [input] = obj.buildVarInput(varargin);
        end
        
    end
    
    methods (Access=private)
        function [input] = buildVarInput(obj,varargin)
            N = numel(varargin);
            base = 'x';
            input = '(';
            for e = 1:N
                input = [input base num2str(e) ','];
            end
            input(end) = ')';
        end
    end
    
end