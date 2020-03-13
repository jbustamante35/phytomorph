classdef storageBlock <  matlab.mixin.Copyable & matlab.mixin.Heterogeneous & handle

    properties

        data;

    end


    methods

        function [obj] = storageBlock(data)
            if nargin == 0
                data = [];
            end
            obj.data = data;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % as if the storageBlock contains input data
        % recursive definition
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [b] = containsInputData(obj)


            if obj.isInputData()
                b = true;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if the data is not all inputBlock type
            % look at each element
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif isa(obj.data,'storageBlock')


                if ~obj.data.isInputData()
                    % loop over the height of the matrix
                    for h = 1:size(obj.data,1)
                        % loop over the width of the matrix
                        for w = 1:size(obj.data,2)
                            % if the [h,w] element is inputData return true
                            if isa(obj.data(h,w),'inputBlock')
                                b(h,w) = true;
                            % if the [h,w] element is storageBlock check
                            elseif isa(obj.data(h,w),'storageBlock')
                                b(h,w) = obj.data(h,w).containsInputData();
                            % is not the above - fail
                            else
                                b(h,w) = false;
                            end
                        end
                    end

                    b = all(b(:));

        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if the data is all inputBlock type
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    b = true;
                end
            


            else
                b = false;


            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        function [b] = isInputData(obj)
            b = isa(obj,'inputBlock');
        end


    end

    methods (Access = protected)

        function [cpObj] = copyElement(obj)
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            if isa(cpObj.data,'storageBlock')
                cpObj.data = copy(cpObj.data);
            end
        end


    end


    methods (Static)
        function [newData] = repmatWide(dataBlock,width)
            for h = 1:size(dataBlock,1)
                for w = 1:width
                    newData(h,w) = copy(dataBlock(h));
                end
            end
        end
    end

end

