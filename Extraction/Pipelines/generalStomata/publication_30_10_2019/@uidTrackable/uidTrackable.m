classdef uidTrackable < handle

    properties
        name;
        createDate;
        uid;
    end

    methods
        
        function [obj] = uidTrackable(name)
            obj.name = name;
            obj.createDate = datestr(datetime,'FFF_ss_hh_dd_mm_yyyy');
            obj.uid = DataHash([obj.createDate obj.createDate]);
        end
    
        function [nameType] = getUidType(obj)
            nameType = [obj.uid '_' class(obj)];
        end

        %%%%%%%%%%%%%%%%%%%%
        % logger function for logging flow(s)
        %%%%%%%%%%%%%%%%%%%%
        function [] = logger(obj,functionCall)
            %{
            % report flow
            if isa(obj,'ceLayer')
                ceType = obj.ceType;
                msg = ['[' obj.name '@' obj.uid ']'];
                switch ceType
                    case 'bottom'
                         timingBlock('stop');
                    case 'top'
                         timingBlock('l-start',msg);
                end
            else
                 msg = ['[' obj.name '@' obj.uid ']'];
            end
            %}
            %fprintf(['Calling ' functionCall ' on '  class(obj) ':' obj.name '-' obj.uid '\n']);
           
        end
    


    end

end