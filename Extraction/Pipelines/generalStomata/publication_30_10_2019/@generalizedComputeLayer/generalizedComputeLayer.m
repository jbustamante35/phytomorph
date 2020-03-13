classdef (Abstract) generalizedComputeLayer ...
            <  uidTrackable & flowable & computable

    properties
       
    end

    methods
        
        function [obj] = generalizedComputeLayer(name)
            obj@uidTrackable(name);
        end
        
    end

    %{
    methods (Static, Sealed, Access = protected)
        function default_object = getDefaultScalarElement
            default_object = generalizedDeviceLayer('default');
        end
    end
    %}
    methods

        %{
        %%%%%%%%%%%%%%%%%%%%
        % logger function for logging flow(s)
        %%%%%%%%%%%%%%%%%%%%
        function [] = logger(obj,functionCall)

            % report flow
            fprintf(['Calling ' functionCall ' on '  class(obj) ':' obj.name '-' obj.uid '\n']);
            % link current layer to previous layer
            %obj.linkCurrentLayer();
            % set current layer to this
            %obj.setAsCurrentLayer();
        end
        %}

        function [] = setAsCurrentLayer(obj)
            mksqlite(['UPDATE currentLayer SET uid =''' obj.uid ''';']);
        end

        function [] = linkCurrentLayer(obj)
            % select the current layer [fid] for attaching the parent
            p = mksqlite(['SELECT uid FROM currentLayer']);
            try
                mksqlite('INSERT INTO layerTable (uid,pid) VALUES (?,?)',obj.uid,p.uid);
            catch
                %fprintf(['[' obj.uid ',' p.uid '] is already logged.\n']);
            end
        end



      

    end

end