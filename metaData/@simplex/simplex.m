classdef simplex
    properties
        key;
        createTime;
        t;
        ptData;
        atData;
        toRender;
    end

    methods
        function [obj] = simplex(ptData,atData,isBasic,toRender,termValue)
            if nargin == 4
                obj.t = term(isBasic,toRender);
            else
                obj.t = termValue;
            end
            obj.toRender = toRender;


            [K,T] = term.keyTime();
            obj.key = K;

            obj.createTime = T;
            
            obj.ptData = ptData;
            obj.atData = atData;

            if obj.toRender
                mksqlite('INSERT INTO myDataTable VALUES (?,?,?,?,?,?)', {'',K,obj.t.key,T,ptData,atData});
            end
        end

        function [ret] = plus(a,b)
            ptData = cat(1,a.ptData,b.ptData);
            atData = cat(1,a.atData,b.atData);
            newTerm = a.t + b.t;
            ret = simplex(ptData,atData,false,false,newTerm);
        end

        function [ret] = minus(a,b)
            ptData = cat(1,a.ptData,b.ptData);
            atData = cat(1,a.atData,b.atData);
            newTerm = a.t - b.t;
            ret = simplex(ptData,atData,false,false,newTerm);
        end

        function [ret] = mtimes(a,b)
            ptData = cat(1,a.ptData,b.ptData);
            atData = cat(1,a.atData,b.atData);
            newTerm = a.t * b.t;
            ret = simplex(ptData,atData,false,false,newTerm);
        end

    end

end