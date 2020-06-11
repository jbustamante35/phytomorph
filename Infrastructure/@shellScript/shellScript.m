classdef shellScript < handle
    
    properties
        oFile;
        cmdSequence;
    end
    
    methods
        function [this] =  shellScript(oFile)
            if nargin == 0;oFile = '';end
            this.cmdSequence = shellCMD('#!/bin/bash');
            this.oFile = oFile;
        end
        
        function [] = addCommand(this,cmd,n)
            if nargin == 2;n = numel(this.cmdSequence)+1;end
            this.cmdSequence(n) = cmd;
        end
        
        function [] = projectScript(this,oFile)
            fileID = fopen(this.oFile,'w');
            for e = 1:numel(this.cmdSequence)
                this.cmdSequence(e).projectScript(fileID);
            end
        end
    end
end


%{
    testScript = shellScript('~/test.sh');
%}