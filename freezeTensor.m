function [fT] = freezeTensor(T,isSet)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this will "freeze" a set of tensors
        % the direct sum analogous to the set
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init return
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fT = [];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % place tensor in cell
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~iscell(T)
            T = {T};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measure grade - the grade from a graded algeb
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        grade = numel(T);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store set size
        % the tensor is [s1,s2...sn|S] or [s1,s2...sn]
        % note that the grade check will "enforce" this
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nargin == 1 || isSet == false
            isSet = false;
            setSize = 1;
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop over each cell/grade to extract the "set" dim
            for g = 1:grade
                tmpSZ = size(T{g});
                grade_check(g) = tmpSZ(end);
            end
            if ~all(grade_check == grade_check(1))
                return 
            end
            setSize = grade_check(1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % handle type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for g = 1:grade
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if isSet then the last index of the array/cell is a "name" 
            % variable
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isSet
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if the grade g is a cell then assume strings else an
                % array of double
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if iscell(T{g})

                    %%%%%%%%%%%%%%%%%%%%
                    % the first object in the cell is 
                    % the prototypical type
                    type{g} = class(T{g}{1});



                    %%%%%%%%%%%%%%%%%%%%
                    if strcmp(type{g},'char')


                        rnk(g) = 2;

                        %%%%%%%%%%%%%%%%%%%%
                        % loop over each element in the 
                        % grade and report the numel from double
                        N = [];
                        for e = 1:numel(T{g})
                            N(e) = numel(double(T{g}{e}));
                        end
                        % max of N - longeset string - not the best solution
                        n = max(N);
                        % dm{g}

                        %%%%%%%%%%%%%%%%%%%%
                        % convert each string to double array of size n
                        Z = zeros(n,numel(T{g}));
                        for e = 1:numel(T{g})
                            tmp = double(T{g}{e});
                            Z(1:numel(tmp),e) = tmp(:);
                        end
                        T{g} = Z;




                        % assume that the grade g is an array of doubles
                        %%%%%%%%%%%%%%%%%%%%
                        % rank of tensor
                        rnk(g) = ndims(T{g});
                        % handle last rank is set index
                        rnk(g) = rnk(g) - 1;
                        % get the dimensions of the tensor
                        dm{g} = size(T{g});
                        % because its a set
                        dm{g}(end) = [];


                        % get pre-shape size
                        tmpSZ = size(T{g});

                        % make the header
                        header{g} = repmat([rnk(g);dm{g}(:);encodeType(type{g})],[1 tmpSZ(end)]);

                        % if set then reshape different
                        newSZ = [numel(T{g})/tmpSZ(end) tmpSZ(end)];
                        % perform reshape
                        T{g} = reshape(T{g},newSZ);




                         %%%%%%%%%%%%%%%%%%%%




                    else
                        fprintf('Problem');
                    end
                    %%%%%%%%%%%%%%%%%%%%



                else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % assume that we have an array of doubles - YES SET
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%
                    % assume that the grade g is an array of doubles
                    %%%%%%%%%%%%%%%%%%%%
                    % rank of tensor
                    rnk(g) = ndims(T{g});
                    % handle last rank is set index
                    rnk(g) = rnk(g) - 1;
                    % get the dimensions of the tensor
                    dm{g} = size(T{g});
                    %because its a set
                    dm{g}(end) = [];

                    % get pre-shape size
                    tmpSZ = size(T{g});

                    % make the header
                    header{g} = repmat([rnk(g);dm{g}(:);encodeType(class(T{g}))],[1 tmpSZ(end)]);

                    % if set then reshape different
                    newSZ = [numel(T{g})/tmpSZ(end) tmpSZ(end)];
                    % perform reshape
                    T{g} = reshape(T{g},newSZ);
                end



            else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % only a single name/number variable NOT A SET
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get the type
                type{g} = class(T{g});
                if strcmp(type{g},'char')
                    % perform conversion
                    T{g} = double(T{g})';
                    % set the rank to one
                    rnk(g) = 2; % odd that all is two
                    % get the size after conversion
                    dm{g} = size(T{g});
                else
                    % if not a string - get the rank
                    rnk(g) = ndims(T{g});
                    % get the size of the tensor 
                    dm{g} = size(T{g});
                end
                % coded - one double 0=double 1=string add more later
                % build the grades header
                header{g} = [rnk(g);dm{g}(:);encodeType(type{g})];
                % call reshape
                T{g} = reshape(T{g},[numel(T{g}) 1]);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store grade
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fT = repmat(grade,[1 setSize]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for each grade
        for g = 1:grade
            % rank of tensor
            %rnk(g) = ndims(T{g});

            % handle last rank is set index
            %if isSet;rnk(g) = rnk(g) - 1;end

            % dim of tensor
            %dm{g} = size(T{g});
            %if isSet;dm{g}(end) = [];end
            % tmp size
            %tmpSZ = size(T{g});
            %newSZ = prod(tmpSZ);
            % if set then reshape different
            %if isSet;newSZ = [newSZ/tmpSZ(end) tmpSZ(end)];end
            %newSZ = [newSZ 1];
            % coded - one double 0=double 1=string add more later
            % call reshape
            %tmpV = reshape(T{g},newSZ);
            % make header - repmat
            % adding 31.10.2019 - type to header
            % header format - [rank [d1 d2..d_rnk]];
            %header = repmat([rnk(g);dm{g}(:)],[1 size(tmpV,2)]);
            % linearize
            fT = [fT;header{g};T{g}];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch ME
        ME
    end
end



function [t] = encodeType(type)
    if strcmp(type,'double')
        t = 0;
    elseif strcmp(type,'char')
        t = 1;
    elseif strcmp(type,'uint16')
        t = 0;
    end
end

%{
T = {rand(4,5),rand(3,4,5,6)};
fT = freezeTensor(T);

T = {rand(1,2),'hello',rand(3,4,5,6)};
fT = freezeTensor(T);
tT = thawTensor(fT);



%}













