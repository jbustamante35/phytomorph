classdef t < handle & matlab.mixin.Heterogeneous
    
    properties
        % store the order and the size
        order = {};
        sz = {};
        signature;
        
        % the tensor to wrap
        T;
        % log of commands
        L = {};
        % log of undo commands
        uL = {};
        
        
        
        % basis vector logs
        bv = {};
        % mean logs
        u = {};
        
        
        
    end
    
    methods
    
        function [obj] = t(T)
            
            obj.T = T;
            
            obj.setSize();
            obj.setOrder();
            
            %{
            sig = ',';
            for e = 1:obj.order{end}
                sig = [sig 'x' num2str(e) ','];
            end
            obj.signature{1} = sig;
            %}
            
            
            %{
            for e = 1:obj.order{end}
                sig{e} = ['x' num2str(e)];
            end
            obj.signature{end+1} = sig;
            %}
        end
        
        function [] = setSize(obj)
            obj.sz{end+1} = size(obj.T);
        end
        
        function [] = setOrder(obj)
            obj.order{end+1} = numel(obj.sz{end});
        end
        
        function [] = updateSignature(obj,sigV)
            obj.signature{end+1} = obj.signature{end}(sigV);
        end
        
        function [] = op(obj,cmd)
            str = strfind(cmd,'{');
            stp = strfind(cmd,'}');
            cmd(1) = [];
            cmd(end) = [];
            cmd = [',' cmd ','];
           
            if obj.toOp(cmd)
                cmd = obj.convertCmd(cmd);
                obj.pr(cmd);
                %obj.p(cmd{1});
                %obj.r(cmd{2});
            end
        end
        
        function [] = func(obj,func)
            for e = 1:numel(obj)
                obj(e).T = func(obj(e).T);
            end
        end
        
        
        
        function [b] = toOp(obj,newOrder)
            b = true;
          %  b = ~strcmp(obj.signature{end},newOrder);
        end
        
        function [b] = doesOrderMatch(obj,cmd)
            b = true;
            %b = numel(strfind(cmd,'x')) == obj.order{end};
        end
        
        
        
        
        
        function [retCmd] = convertCmd(obj,cmd)
           
            newCmd = [];
            pCmd = [];
            rCmd = [];
            
            ridx = strfind(cmd,'r');
            redCmd = cmd(ridx+1);
            argStr = strfind(cmd,'(');
            argStp = strfind(cmd,')');
            reduceCmd = [];
            for e = 1:numel(ridx)
                tmpIdx1 = ridx(e) < argStr;
                tmpIdx2 = argStr(tmpIdx1(1)) < argStp;
                tmpIdx1 = argStr(tmpIdx1(1));
                tmpIdx2 = argStp(tmpIdx2(1));
                reduceCmd(1,e) = str2num(redCmd(e));
                reduceCmd(2,e) = str2num(cmd((tmpIdx1+1):(tmpIdx2-1)));
            end
            
            for e = 1:numel(argStr)
                cmd(argStr(e):argStp(e)) = [];
                argStr = strfind(cmd,'(');
                argStp = strfind(cmd,')');
            end
            
            
            cmd = strrep(cmd,'r','x');
            
            
            
            
            cidx = strfind(cmd,',');
            for e = 1:(numel(cidx)-1)
                str = cidx(e)+1;
                stp = cidx(e+1)-1;
                cmdSnip = cmd(str:stp);
                
                
                
                
                cmdSnip = [cmdSnip 'x'];
                xidx = strfind(cmdSnip,'x');
                miniCmd = [];
                for x = 1:(numel(xidx)-1)
                    subSnip = cmdSnip((xidx(x)+1):(xidx(x+1)-1));
                    pCmd = [pCmd str2num(strrep(subSnip,'x',''))];
                    miniCmd = [miniCmd str2num(strrep(subSnip,'x',''))];
                end
                
                %rCmd(e) = prod(obj.sz{end}(miniCmd));
                rCmd{e} = miniCmd;
                
                
                %cmdSnip = strrep(cmdSnip,'x','');
                
                %newCmd = [newCmd str2num(cmdSnip)];
            end
            
            here = 1;
            retCmd{1} = reduceCmd;
            retCmd{2} = pCmd;
            retCmd{3} = rCmd;
        end
        
        
        % permute command
        function [] = p(obj,cmd)
            % generate commands
            fCommand = obj.genDo('permute',cmd);
            rCommand = obj.genunDo('permute',cmd);
            % execute comamnds
            obj.T = fCommand(obj.T);
            % log commands
            obj.L{end+1} = {fCommand};
            obj.uL{end+1} = {rCommand};
            % set order and size
            obj.setSize();
            obj.setOrder();
            %obj.updateSignature(cmd);
        end
        
        function [] = r(obj,cmd)
            % generate commands
            fCommand = obj.genDo('reshape',cmd);
            rCommand = obj.genunDo('reshape',cmd);
            % execute comamnds
            obj.T = fCommand(obj.T);
            % log commands
            obj.L{end+1} = {fCommand};
            obj.uL{end+1} = {rCommand};
            % set order and size
            obj.setSize();
            obj.setOrder();
            %obj.updateSignature(cmd);
        end
        
        function [] = pr(obj,cmd)
            
            % generate forward command for pca
            red_fCommand = obj.genDo('pca',cmd{1});
            
            % run pca on needed dims
            [out] = red_fCommand(obj.T);
            obj.bv{end+1} = out.bv;
            obj.u{end+1} = out.u;
            obj.T = out.data;
            
            % generate the reverse pca
            red_rCommand = obj.genunDo('pca',cmd{1});
           
             % set order and size
            obj.setSize();
            obj.setOrder();
            
            % generate the rehape co
            sz = size(obj.T);
            reshapeArgs = cmd{3};
            for e = 1:numel(reshapeArgs)
                narg(e) = prod(sz(reshapeArgs{e}));
            end
            
            % generate commands
            per_fCommand = obj.genDo('permute',cmd{2});
            per_rCommand = obj.genunDo('permute',cmd{2});
            
            % execute comamnds
            obj.T = per_fCommand(obj.T);
            
            
            % set order and size
            obj.setSize();
            obj.setOrder();
            
            
            
            
            % generate commands
            re_fCommand = obj.genDo('reshape',narg);
            % execute comamnds
            obj.T = re_fCommand(obj.T);
            re_rCommand = obj.genunDo('reshape',cmd{3});
            
            
            % set order and size
            obj.setSize();
            obj.setOrder();
           
            % log commands
            obj.L{end+1} = {red_fCommand per_rCommand re_fCommand};
            obj.uL{end+1} = {re_rCommand per_rCommand red_rCommand};
        end
        
       
        
        
        
        
        function [cmd] = genDo(obj,type,args)
            switch type
                case 'permute'
                    cmd = @(X)permute(X,args);
                case 'reshape'
                    cmd = @(X)reshape(X,args);
                case 'pca'
                    cmd = @(X)t.d(X,args);
            end
        end
        
        function [cmd] = genunDo(obj,type,args)
            switch type
                case 'permute'
                    cmd = @(X)ipermute(X,args);
                case 'reshape'
                    values = obj.sz{end};
                    cmd = @(X)reshape(X,values);
                case 'pca'
                    values.bv = obj.bv{end};
                    values.u = obj.u{end};
                    cmd = @(X)t.ud(X,args,values);
            end
        end
        
        function [] = undo(obj)
            cmdSequence = obj.uL{end};
            for c = 1:numel(cmdSequence)
                obj.T = cmdSequence{c}(obj.T);
                obj.cleanLastProperties();
            end
            obj.cleanLastCommands();
           
        end
        
        function [] = cleanLastProperties(obj)
            obj.sz(end) = [];
            obj.order(end) = [];
            %obj.signature(end) = [];
        end
        
        function [] = cleanLastCommands(obj)
            obj.L(end) = [];
            obj.uL(end) = [];
        end
        
        
        %{
        function [varargout] = subsref(obj,S)
            
            if strcmp(S(1).type,'()')
                if nargout >= 1
                    [varargout{1:nargout}] = builtin('subsref',obj.T,S);
                else
                    builtin('subsref',obj.T,S)
                end
            else
                if nargout >= 1
                    [varargout{1:nargout}] = builtin('subsref',obj,S);
                else
                    builtin('subsref',obj,S)
                end
            end
              
        end
        %}
        
    end
    
    
    methods (Static)
   
        function [out] = d(X,cmd)
            if ~isempty(cmd)
                dimOp = cmd(1,:);
                dimRe = cmd(2,:);
                bv = {};
                u = {};
                for e = 1:numel(dimOp)
                    pvec = [dimOp(e) setdiff(1:ndims(X),dimOp(e))];
                    tmp = permute(X,pvec);
                    tsz = size(tmp);
                    tmp = reshape(tmp,[tsz(1) prod(tsz(2:end))]);
                    [U,E,L] = PCA_FIT_FULL_Tws(tmp,dimRe(e));
                    S = PCA_REPROJ_T(tmp,E,U);
                    tsz(1) = dimRe(e);
                    tmp = reshape(S,tsz);
                    tmp = ipermute(tmp,pvec);
                    X = tmp;
                    bv{end+1} = E;
                    u{end+1} = U;
                end
                out.data = X;
                out.bv = bv;
                out.u = u;
            else
                out.data = X;
                out.bv = 1;
                out.u = 0;
            end
         end
        
         
        function [X] = ud(X,cmd,args)
            if ~isempty(cmd)
                dimOp = cmd(1,:);
                dimRe = cmd(2,:);
                for e = 1:numel(dimOp)
                    pvec = [dimOp(e) setdiff(1:ndims(X),dimOp(e))];
                    tmp = permute(X,pvec);
                    tsz = size(tmp);
                    tmp = reshape(tmp,[tsz(1) prod(tsz(2:end))]);
                    X = PCA_BKPROJ_T(X,args.bv{e},args.u{e});
                    tsz(1) = size(X,1);
                    X = reshape(X,tsz);
                    X = ipermute(X,pvec);
                end
            else
                X = X;
            end
        end
        
         
    end
    
end

%{


    T = rand(30,4,100);
    g = t(T);
    g.op('{x2x1,x3}');

 

    T = rand(30,4,100);
    g = t(T);
    g.op('{r1(2),x3,x2}');


    r = rand(1,4);
    g = t(r);
    g = repmat(g,[100 4]);
    g.op('{x3,r1(2),x2}');



    T = rand(30,4,100);
    g = t(T);
    g.op('{x3,r1(2),x2}');


    T1 = rand(30,4,100);
    T2 = rand(30,5,100);
    g1 = t(T1);
    g2 = t(T2);
    G = [g1 g2];

    func = @(X)squeeze(sum(X,2));
    G.func(func);

%}

