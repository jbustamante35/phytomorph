classdef addressBookEntry < doid

    properties
        mainName;
        pseudonymList;
        channel;
    end
    
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%
        % entry in book
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = addressBookEntry(mainName,channel,pseudonymList)
            if nargin ~= 0
                obj.mainName = mainName;
                obj.channel = channel;
            end
            % pseudonym list is a cell array of strings
            % monikers that the entry "goes" by
            if nargin == 3
                obj.pseudonymList = pseudonymList;
            else
                obj.pseudonymList = {};
            end
          
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        % add to the list of names
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [] = addpseudonym(obj,pseudonym)
            obj.pseudonymList{end+1} = pseudonym;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        % is a msg to this entry
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [b] = isTo(obj,name)
            b = strcmp(name,'*');
            b = b | contains(obj.pseudonymList,name);
            b = b | strcmp(name,obj.mainName);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % check for equality
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [b] = eq(obj,a)
            b = strcmp(obj.mainName,a.mainName);
        end
        
        function [] = combine(obj,a)
            if obj == a
                obj.pseudonymList = unique({obj.pseudonymList{:},a.pseudonymList{:}});
            end
        end
    end
    
end