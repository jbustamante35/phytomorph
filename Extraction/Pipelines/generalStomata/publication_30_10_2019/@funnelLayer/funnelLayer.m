classdef funnelLayer < generalizedFunctionLayer

    properties
        commandWidth;
        glueCommands;
    end


    methods
    
        function [obj] = funnelLayer(name,glueCommands)
            obj = obj@generalizedFunctionLayer(name);
            obj.commandWidth = numel(glueCommands);
            obj.glueCommands = glueCommands;
        end

        function [result] = compute(obj,varargin)
            result = varargin;
        end

        function [] = configure(obj)
        end

    end

    methods (Access=private)


        function [out] = glue(obj,in1,in2)
            for e = 1:numel(in1)
                curType = obj.glueCommands{e};
                if strcmp(curType{1},'cat')
                    out{e} = catType(in1{e},in2{e},curType{2});
                elseif strcmp(curType{1},'keep')
                    out{e} = catType(in1{e},in2{e},curType{2});
                end
            end
        end


        function [out] = catType(in1,in2,dim)
            out= cat(dim,in1,in2);
        end

        function [out] = keepType(in1,in2,select)
            if (select == 1)
                out = in1;
            elseif (select == 2)
                out = in2;
            end
        end

    end


    methods (Static)
        
        function [cmd] = defaultCommands(width)
            for e = 1:width
                cmd{e} = {'cat',1};
            end
        end

    end

end