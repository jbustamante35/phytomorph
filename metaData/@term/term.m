classdef term < handle
    properties

        data = '';
        key = '';
        createTime = '';
        toRender;
    end
    methods
        function [obj] = term(isBasic,toRender)


            obj.toRender = toRender;


            [K,T] = term.keyTime();
            obj.key = K;

            obj.createTime = T;

            if isBasic
                obj.data = ['|' K '|'];
            else
                obj.data = '';
            end

            if obj.toRender 
                data = [obj.data];
                type = 'term';
                mksqlite('INSERT INTO myStringTable VALUES (?,?,?,?,?)', {'',K,T,data,type});
            end

        end

        function [n] = numelTerms(obj)
            strRep = strrep(obj.data,'|','+');
            fidx = [strfind(strRep,'+') strfind(strRep,'-')];
            n = numel(fidx) - 1;
        end

        function [n] = numelNodes(obj)
            strRep = strrep(obj.data,'|','*');
            fidx = [strfind(strRep,'*')];
            n = numel(fidx) - 1;
        end

        function [ret] = removeNode(obj,n)
            [tmpData] = term.removeTerm(obj.data,n);
            ret = term(false,false);
            ret.data = tmpData;
        end

        function [ret,sgn] = getTerm(obj,n)
            strRep = obj.data;
            [ret,sgn] = term.getTermFromExpression(strRep,n);
        end

        function [res] = d(obj)
            N = obj.numelTerms();
            res = term(false,false);
            for termN = 1:N
                restmp = term(false,false);
                [tmp,sgn] = obj.getTerm(termN);
                tmpTerm = term(false,false);
                tmpTerm.data = tmp;
                for node = 1:tmpTerm.numelNodes()
                    if rem(node,2) == sgn
                        restmp = restmp + tmpTerm.removeNode(node);
                    else
                        restmp = restmp - tmpTerm.removeNode(node);
                    end
                end
                res = res + restmp;
            end
        end


        function [ret] = mtimes(obj,b)
            aStr = obj.toString;
            bStr = b.toString;
            ret = term(false,false);
            ret.data = ['|' aStr '*' bStr '|'];
        end


        function [ret] = plus(a,b)
            aStr = a.toString;
            if ~isempty(aStr)
                aStr = [aStr '+'];
            end
            bStr = b.toString;
            if isempty(bStr)
                if ~isempty(aStr)
                    aStr(end) = [];
                end
            else
                if strcmp(bStr(1),'-')
                    aStr(end) = [];
                end
            end
            ret = term(false,false);
            ret.data = ['|' aStr bStr '|'];
        end


        function [ret] = minus(a,b)
            aStr = a.toString;
            bStr = b.toString;
            if ~isempty(bStr)
                bStr = ['-' bStr];
            end
            ret = term(false,false);
            ret.data = ['|' aStr bStr '|'];
        end

        function [str] = toString(obj)
            str = obj.data;
            if ~isempty(str)
                str(1) = [];
                if ~isempty(str)
                    str(end) = [];
                end
            end
%{
            str = '';
            for e = 1:numel(obj)
                tmp = obj(e).data;
                tmp(1) = [];
                tmp(end) = [];
                str = [str '|' tmp];
            end
            str = [str '|'];
%}
        end
        function [] = render(obj)
            type = 'term';
            sql = ['UPDATE myStringTable SET data = ''' obj.data ''' WHERE key=''' obj.key ''''];
            q = mksqlite(sql);
        end


        
    end

    methods (Static)
        function [] = query()
            sql = ['SELECT * FROM myStringTable'];
            q = mksqlite(sql);
            for e = 1:numel(q)
                q(e)
            end
        end

        function [ret] = removeTerm(strRep,n)
            strRep = strrep(strRep,'|','*');
            glueIDX = strfind(strRep,'*');
            strRep((glueIDX(n)+1):glueIDX(n+1)) = [];
            ret = strRep;
            ret(1) = '|';
            ret(end) = '|';
        end


        function [ret,sgn] = getTermFromExpression(strRep,n)
            ogn = strRep;
            strRep = strrep(strRep,'|','+');
            strRep = strrep(strRep,'-','+');
            glueIDX = strfind(strRep,'+');
            ret = strRep((glueIDX(n)):(glueIDX(n+1)));
            ret = strrep(ret,'+','|');
            sgn = ogn(glueIDX(n));
            sgn = ~strcmp(sgn,'-');
        end
        
        function [K,T] = keyTime()
            % generate random index key
            K = randi((2^32)-1,1,'uint32');
            K = num2str(K);
            T = datenum(datetime);
        end

        function [] = closeDB()
            mksqlite('close');
        end

        function [] = openDB()
            mksqlite('open',':memory:');
            % set blob type
            mksqlite('typedBLOBs', 1 );
            % wrapping
            mksqlite('param_wrapping', 1 );
            mksqlite('CREATE TABLE myStringTable (id INTEGER PRIMARY KEY AUTOINCREMENT,key,time,data,type)');
            mksqlite('CREATE TABLE myDataTable (id INTEGER PRIMARY KEY AUTOINCREMENT,key,termKey,time,ptData,atData)');
        end


    end

end