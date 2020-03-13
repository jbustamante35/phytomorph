classdef algorithm < layerStack & testable

    properties
    
        algorithmName;
        author;
        algorithmVersion;
        description;
        contributorList;
        lastTestDate;
    end
    
    methods 
    
        function [obj] = algorithm(name,layers)
            obj@layerStack(name,layers);
            obj@testable(name);
            obj.topCE.bannerMessage{1} = ['Algorithm Name: ' name];
        end
        
        function [] = setVersion(obj,version)
            obj.algorithmVersion = version;
            obj.topCE.bannerMessage{2} = ['Version: ' version];
        end
        
        function [] = setAuthor(obj,author)
           obj.author = author;
           obj.topCE.bannerMessage{3} = ['Author: ' author];
        end
        
        function [] = setDescription(obj,description)
           obj.description = description;
           obj.topCE.bannerMessage{4} = ['Description: ' description];
        end
        
        function [] = attachNote(obj,contributor,noteString)
            obj.contributorList{end+1} = contributor;
            obj.contributorList = unique(obj.contributorList);
        end
        
        function [] = incrementVersion(obj)
            obj.algorithmVersion = obj.algorithmVersion + 1;
        end
        
        function [] = setLastTestDate(obj,testDate)
        
        end
        
    end
    
end