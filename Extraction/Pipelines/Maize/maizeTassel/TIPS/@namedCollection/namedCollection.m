classdef namedCollection < doidmm

    properties
        % name of the collection
        name;
        % order of the collection - number of index
        order;
        % the collection container
        data;
    end
    
    methods
        
        function [this] = namedCollection(name,initData,order)
            switch class(initData)
                case 'char'
                    switch initData
                        case 'array'
                            this.data = [];
                        case 'cell'
                            this.data = {};
                        case 'map'
                            this.data = containers.Map;
                        case 'htcResource'
                            this.data = htcResource.empty(0,1);
                        case 'htcComputeNode'
                            this.data = htcComputeNode.empty(0,1);
                        otherwise
                            this.data = initData;
                    end
                otherwise 
                    this.data = initData;
            end
            this.name = name;
            this.order = 1;
            if nargin == 3;this.order = order;end
        end
        
        function [n] = end(this,k,n)
            if k <= this.order
                n = builtin('end',this.data,k,n);
            else
                syms n
                %n=0;
                %n = '1';
            end
            
        end
        
        function [element] = getElement(this,index)
            if isa(this.data,'cell')
                element = this.data{index};
            else
                element = this.data(index,:);
            end
        end
        
        function [varargout] = subsref(this,subs)
            
            %[subs,normalSubs] = buildSubElement(this,subs);
            switch subs(1).type
                case '.'
                    if strcmp(subs(1).subs,'ref')
                        varargout{1:nargout} = this;
                    else
                        [varargout{1:nargout}] = builtin('subsref',this,subs);
                    end
                case '()'
                    try
                        [varargout{1:nargout}] = builtin('subsref',this.data,subs);
                    catch
                        varargout{1} = '';
                    end
                case '{}'
                    [varargout{1:nargout}] = builtin('subsref',this.data,subs);
            end
            %{
            if ~isempty(normalSubs.subs)
                a = subsref(varargout{1},normalSubs);
                if nargout == 0
                    [varargout{1}] = a;
                else
                    [varargout{1:nargout}] = a;
                end
            end
            %}
            here = 1;
        end
        
        function [varargout] = subsasgn(this,subs,in)
             switch subs(1).type
                case '.'
                    [varargout{1:nargout}] = builtin('subsasgn',this,subs,in);
                case '()'
                    [this.data] = builtin('subsasgn',this.data,subs,in);
                case '{}'
                    [this.data] = builtin('subsasgn',this.data,subs,in);
             end
             varargout = {this};
        end
        
        function [subs,normalSubs] = buildSubElement(this,subs)
            normalSubs = subs;
            normalSubs.subs = normalSubs.subs((this.order+1):end);
            
            if ~isempty(normalSubs.subs)
                if isa(normalSubs.subs{1},'cell')
                    normalSubs.type = '{}';
                    here = 1;
                    normalSubs.subs{1} = normalSubs.subs{1}{1};
                else
                    normalSubs.type = '()';
                end
            end
            
            if isa(subs.subs,'cell')
                subs.subs = subs.subs(1:this.order);
            end
        end
        
        function [] = gather(this,func)
            here = 1;
        end
        
        
    end
end