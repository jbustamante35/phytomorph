classdef m95 < doidmm
    
    properties
        % name of function
        name;
        % function to call
        baseFunc; 
        % not sure i need this yet
        signature;
        % arguments for this function
        X;
    end
    
    methods
        
        function [this] = m95(baseFunc,name)
            if nargin == 1;name = '';end
            this.name = name;
            if isa(baseFunc,'char')
                this.signature = formalFunc.signatureFromName(baseFunc);
                this.baseFunc = formalFunc.funcFromName(baseFunc);
            elseif isa(baseFunc,'function_handle')
                this.baseFunc = baseFunc;
            end
            this.signature = formalFunc.signatureFromHandle(this.baseFunc);
            this.X = z0.buildArray(this.signature.size(2));
        end
        
        function [varargout] = subsref(this,subs)
            % switch 
            switch subs(1).type
                case '.'
                    [varargout{1:nargout}] = builtin('subsref',this,subs);
                case '()'
                    
                case '{}'
                    for e = 1:numel(this.X)
                        bool(e) = isa(this.X(e),'z0');
                    end
                    fidx0 = find(bool);
                    fidx1 = find(~bool);
                    % make empty array of z0
                   	for e = 1:numel(this.X)
                        newX(e) = z0;
                    end
                    % for each input argument
                    for e = 1:numel(subs(1).subs)
                        if ~isa(subs(1).subs{e},'fa')
                            % fill the needed spots with the arguments
                            newX(fidx0(e)) = z(subs(1).subs{e});
                        end
                    end
                    % assign the base function
                    varargout{1} = m95(this.baseFunc);
                    % assign the arguments
                    varargout{1}.X = newX;
                    
                    
                    % save the output m95 object
                    varargout{1}.persist();
            end
        end
        
        
        function [] = persist(this)
            here = 1;
        end
        
        
    end
    
    
    methods (Static)
        
        
        
    end
    
    
end