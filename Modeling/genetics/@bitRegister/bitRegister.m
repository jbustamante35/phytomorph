classdef bitRegister
    properties
        data;
        MAX_B;
    end


    methods
        function [obj] = bitRegister(MAX_B,data)
            obj.MAX_B = MAX_B;
            obj.data = data;
        end

        function [ret] = get(obj,n,type)
            ret = obj.data(:,:,n);
        end

        function [seqLength] = getSeqLength(obj)
            seqLength = size(obj.data,3);
        end


        function [X] = entangleRegisterAtN(obj,n,type)
            if nargin == 2
                type = false;
            end
            tmpData = obj.get(n);
            X = bitRegister.entangle(tmpData,obj.MAX_B,type);
        end

        function [entangled] = entangleNewBits(obj,n,newData)
            [ret] = obj.get(n,'dirac');
            ret = [ret newData];
            ret = squeezeQbits2(ret,{[1:obj.MAX_B],[1:size(newData,2)]},3,'uint8');
            entangled = bitRegister.entangle(ret,obj.MAX_B + size(newData,2),false);
        end
        
        
    end


    methods (Static)
        
        function [X] = entangle(X,MAX_B,type)
            [X] = quantumEncode(MAX_B,X,type);
        end
    end
end