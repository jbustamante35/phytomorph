classdef (Abstract) generalizedComputeLayer < matlab.mixin.Heterogeneous & handle

    properties
        name;
        createDate;
        uid;
        %%%
        parentList;
        childList;
    end

    methods
        
        function [obj] = generalizedFeatureLayer(name)
            obj.name = name;
            obj.createDate = datestr(datetime,'ss_hh_dd_mm_yyyy');
            obj.uid = DataHash([obj.createDate obj.createDate]);
        end
        
    end
    
    methods (Abstract)
        
        feature = computeFeature(obj,image,point,domain,varargin);
        
        
    end
end