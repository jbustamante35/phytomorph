function [nLevel] = timingBlock(type,blockMsg)
    global queue
    breakBlockLength = 150;
    switch type
        case 'l-start'
            nLevel = numel(queue)+1;
            msg = blockMsg;
            spaceBlock = repmat('   ',[1 nLevel-1]);
            r = repmat('---',[1 nLevel-1]);
            lb = repmat('-',[1 30]);
            %fprintf(['|' lb '\n']);
            %fprintf(['|' r '>[START LAYER]-' msg '\n']);
            lb = repmat('-',[1 breakBlockLength]);
            fprintf([spaceBlock '|' lb '\n']);
                 
            fprintf([spaceBlock '|-' '>[START LAYER]-' msg '\n']);
            spaceBlock = repmat('   ',[1 nLevel]);
            fprintf([spaceBlock '|' '\n']);
            %fprintf([spaceBlock '|__' '\n']);
            queue{end+1} = {type,blockMsg,clock};
        case {'b-start','start'}
            nLevel = numel(queue)+1;
            spaceBlock = repmat('   ',[1 nLevel-1]);
            msg = blockMsg;
            r = repmat('---',[1 nLevel]);
            %fprintf(['|' r '>[START]-' msg '\n']);
            fprintf([spaceBlock '|' '-' '>[START BLOCK]-' msg '\n']);
            spaceBlock = repmat('   ',[1 nLevel]);
            fprintf([spaceBlock '|' '\n']);
            %fprintf([spaceBlock '|__' '\n']);
            queue{end+1} = {type,blockMsg,clock};
        case 'stop'
            nLevel = numel(queue);
            spaceBlock = repmat('   ',[1 nLevel-1]);
            
            if nLevel > 0
                type = queue{end}{1};
               
                if nargin == 1;msg = queue{end}{2};end
                if nargin == 2;msg = blockMsg;end
                tm = queue{end}{3};
                deltaTM = etime(clock,tm);
                
                switch type
                    case 'l-start'
                        r = repmat('---',[1 nLevel]);
                        %fprintf(['|' r '>[STOP ]-' msg '-[' num2str(deltaTM) '(s)]\n']);
                        %fprintf([spaceBlock '___|' '\n']);
                        %fprintf([spaceBlock '|' '\n']);
                        EspaceBlock = repmat('   ',[1 nLevel]);
                        fprintf([EspaceBlock '|' '\n']);
                        fprintf([spaceBlock '|' '-' '>[STOP LAYER]-' msg '-[' num2str(deltaTM) '(s)]\n']);
                        
                        lb = repmat('-',[1 breakBlockLength]);
                        fprintf([spaceBlock '|' lb '\n']);
                    case {'b-start','start'}
                        r = repmat('---',[1 nLevel]);
                        %fprintf(['|' r '>[STOP ]-' msg '-[' num2str(deltaTM) '(s)]\n']);
                        %fprintf([spaceBlock '___|' '\n']);
                        %fprintf([spaceBlock '|' '\n']);
                        EspaceBlock = repmat('   ',[1 nLevel]);
                        fprintf([EspaceBlock '|' '\n']);
                        fprintf([spaceBlock '|' '-' '>[STOP]-' msg '-[' num2str(deltaTM) '(s)]\n']);
                        %fprintf([spaceBlock '|' '\n']);
                end        
                queue(end) = [];
                
              
            else
                fprintf(['TIME BLOCK error.\n']);
            end
            
            
        case 'note'
            nLevel = numel(queue);
            spaceBlock = repmat('   ',[1 nLevel]);
           
            if nLevel > 0
                msg = blockMsg;
                tm = queue{end}{3};
                deltaTM = etime(clock,tm);
                r = repmat('   ',[1 nLevel]);
                fprintf([spaceBlock '|-' '>[NOTE]-' msg '-[' num2str(deltaTM) '(s)]\n']);
            else
                fprintf(['TIME BLOCK error.\n']);
            end
            
            
        case 'clear'
            queue = {};
            
        case 'levelQuery'
            nLevel = numel(queue);
    end

end