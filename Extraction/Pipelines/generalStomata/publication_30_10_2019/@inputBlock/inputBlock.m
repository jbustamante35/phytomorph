classdef inputBlock < storageBlock


    properties

    	lastLayer;      % last layer on record

    end

    methods

        function [obj] = inputBlock(data)
            if ~isa(data,'cell')
                fprintf(['Warning as of 14.11.2019: Data presented to inputBlock was not a cell.\nds'])
                %data = {data};
            end
            obj = obj@storageBlock(data);
            obj.lastLayer = 'null';
        end

        function [out] = getInput(obj)
            out = obj.data;
        end


        function [] = setLastLayer(obj,lastLayer)
            obj.lastLayer = lastLayer;
        end


    end

end