classdef shellCMD < doid

    properties
        cmtString;
        cmdString;
    end

    methods
        
        function [this] = shellCMD(cmdString,comment)
            if nargin == 1;comment = '';end
            this.cmdString = cmdString;
            this.cmtString = comment;
        end
        
        function [] = projectScript(this,oFile)
            if nargin == 1;oFile = this.oFile;end
            if isa(oFile,'char');oFile = fopen(oFile,'w');end
            % write the comment string if not empty
            if ~isempty(this.cmdString)
                fprintf(oFile,this.cmdString);
            end
            % write the command
            fprintf(oFile,this.cmdString);
        end
        
        
    end
    
end