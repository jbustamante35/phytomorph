classdef rGrid < handle
    
    properties
        % set of state vectors for -- [folded/unfolded base/apex];
        state;
        size;
        data;
        grid;
        sz;
        vFunc;
    end

    
    methods
        function [obj] = rGrid(minLength,maxLength,npLength,minWidth,maxWidth,npWidth,vFunc)
            % render the grid
            [w,l] = ndgrid(linspace(minWidth,maxWidth,npWidth),linspace(minLength,maxLength,npLength));
            % size of the stack of data
            obj.sz = [size(l) 6];
            % create the time information
            t = zeros(numel(l),1);
            % create the angle informtion
            a = zeros(size(t));
            % attach the velocity
            obj.attachVelocityFunction(vFunc);
            
            % stack the information [length,width,time,angle,X,Y];
            obj.grid = [l(:),w(:),t(:),a(:),l(:),w(:)];
            % set state
            obj.state = {'fold','base'};
            % get init velocity
            vl = obj.getVelocity('base');
            vw = zeros(size(vl));
            
            % get the frame
            F = obj.getFrame();
            % init the v TAN
            vTAN = bsxfun(@times,vl,F(:,:,:,1));
            % init the vNOR
            vNOR = bsxfun(@times,vw,F(:,:,:,2));
            % init the vx
            vx = vTAN(:,:,1) + vNOR(:,:,1);
            % init the vy
            vy = vTAN(:,:,2) + vNOR(:,:,2);
            
            vl = vx;
            vw = vy;
            
            % include the current velocity as a property of the grid
            obj.grid = [obj.grid,vl(:),vw(:),vx(:),vy(:)];
            % reset the size of the information package
            obj.sz = [size(l) 10];
        end
        
        
        function [] = attachVelocityFunction(obj,vFunc)
            obj.vFunc = vFunc;
        end
        
        function [] = updateGrid(obj,n,dt,h,toQuiver)
            LV = [];
            h2 = figure;
            h3 = figure;
            measurement = [];
            for loop = 1:n
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % make measurements
                %what = getMeasurement(obj,'cellProfile');
                %what(loop) = getMeasurement(obj,'tip_angle');
                measurement = obj.getMeasurement('display',measurement,h2);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % render the grid
                obj.render(h,toQuiver);
                
                
                
                
                %%%%
                % get the speed along the midline from the base perspeective
                vl = obj.getVelocity('base');
                % get the speed from the midlien from the apex persepective
                vl_apex = obj.getVelocity();
                % set the width speed to zero
                vw = zeros(size(vl));
                 
                
                
                % unfold the grid
                obj.stateT('unfold');
                % get the length,width,x and y
                l = obj.grid(:,:,1);
                w = obj.grid(:,:,2);
                x = obj.grid(:,:,5);
                y = obj.grid(:,:,6);
                vl = obj.grid(:,:,7);
                vw = obj.grid(:,:,8);
                vx = obj.grid(:,:,9);
                vy = obj.grid(:,:,10);
              
                
                [~,DX] = obj.getFrame();
                
                [dvxdl,dvxdw] = gradient(vx);
                dvxdl = dvxdl.*DX(:,:,1).^-1;
                dvxdw = dvxdl.*DX(:,:,2).^-1;
                
                
                [dvydl,dvydw] = gradient(vy);
                dvydl = dvydl.*DX(:,:,1).^-1;
                dvydw = dvydw.*DX(:,:,2).^-1;
                
                
                
                
                [ds] = gradient(l);
                
               
                % calculate the increase in length along the midline for
                % each point - from the base perspetive
                dl = vl*dt;
                
                
                % perspective of the length to base
                obj.stateT('base');
                % increase the length
                obj.grid(:,:,1) = obj.grid(:,:,1) + dl;
                obj.grid(:,:,1) = obj.calcLengthFromGrid();
                
                
                
                parfor pt = 1:numel(x)
                    mask = l <= l(pt);
                    mask_vx = vx.*mask;
                    mask_vy = vy.*mask;
                    
                    %dVX(pt) = sum(mask_vx(:)) / sum(mask(:));
                    %dVY(pt) = sum(mask_vy(:)) / sum(mask(:));
                    
                    dVX(pt) = sum(mask_vx(:));
                    dVY(pt) = sum(mask_vy(:));
                    
                    
                    if isnan(dVX(pt))
                        dVX(pt) = 0;
                    end
                    if isnan(dVY(pt))
                        dVY(pt) = 0;
                    end
                end
                dVX = reshape(dVX,obj.sz(1:2));
                dVY = reshape(dVY,obj.sz(1:2));
                here = 1;
                
                
                %dVX = dl;
                %dVY = vw;
                
                
                % get the TAN-NOR frame(s) at each point
                %F = obj.getFrame();
                
                
                % multiply the TAN by the speed along the 
                %vx = bsxfun(@times,vl,F(:,:,:,1));
                %vy = bsxfun(@times,vw,F(:,:,:,2));
                %v = vx + vy;
                
                
                
                
                 
                
                %{
                vw = zeros(size(vl));
                F = obj.getFrame();
                vx = bsxfun(@times,vl,F(:,:,:,1));
                vy = bsxfun(@times,vw,F(:,:,:,2));
                vT = vx + vy;
                [cav] = curl(x,y,vT(:,:,1),vT(:,:,2));
                cav = -cav;
                %}
                
                
                
              
                %cav = -cav;
                %cav = cumsum(cav,2);
                %obj.grid(:,:,4) = obj.grid(:,:,4) + cav*dt;
               
                
                %{
              
                [cav] = curl(x,y,vl_apex,vw);
                parfor pt = 1:numel(x)
                    pt;
                    
                    
                    cur_x = x(pt);
                    cur_y = y(pt);
                    
                    
                    delta_x = x - cur_x;
                    delta_y = y - cur_y;
                    mask = l < l(pt);
                    
                    delta = (delta_x.^2 + delta_y.^2).^.5;
                    
                    delta(pt) = 1;
                    delta_x = delta_x .*delta.^-1;
                    delta_y = delta_y .*delta.^-1;
                    delta(pt) = 0;
                    
                    
                    displace_tan = delta - delta.*cos(cav*dt);
                    displace_nor = delta.*sin(cav*dt);
                    
                    
                    delta_tan_vec = cat(3,delta_x,delta_y);
                    delta_nor_vec = cat(3,delta_y,-delta_x);
                    
                    displace_tan = mask.*displace_tan;
                    displace_nor = mask.*displace_nor;
                    
                    
                    aT = bsxfun(@times,displace_tan,delta_tan_vec);
                    aN = bsxfun(@times,displace_nor,delta_nor_vec);
                    aTOT = aT + aN;
                    
                    rX(pt) = sum(sum(aTOT(:,:,1))) / sum(mask(:));
                    rY(pt) = sum(sum(aTOT(:,:,2))) / sum(mask(:));
                    if isnan(rX(pt))
                        rX(pt) = 0;
                    end
                    if isnan(rY(pt))
                        rY(pt) = 0;
                    end
                end
                rX = reshape(rX,obj.sz(1:2));
                rY = reshape(rY,obj.sz(1:2));
                
                
                
                
                %}
                
                
                
                
                
                %maxL = max(max(obj.grid(:,:,1)));
                %np = size(obj.grid,1);
                
                
                % this will re-space the grid
                % but the [lw] co-ordinates will be wrong
                
                %ds = maxL/np;
                
                %[da] = gradient(obj.grid(:,:,4));
                %[ds] = gradient(obj.grid(:,:,1));
                %K = da.*(ds.^-1);
                
                
                
                %{
                initX = obj.grid(:,1,5);
                initY = obj.grid(:,1,6);
                dx = ds.*cos(-obj.grid(:,:,4));
                dy = ds.*sin(-obj.grid(:,:,4));
                newX = bsxfun(@plus,cumsum(dx,2),initX);
                newY = bsxfun(@plus,cumsum(dy,2),initY);
                obj.grid(:,:,5) = newX;
                obj.grid(:,:,6) = newY;
                %}
                
                
                %update position
                newX = x + dVX*dt;
                newY = y + dVY*dt;
                obj.grid(:,:,5) = newX;
                obj.grid(:,:,6) = newY;
                
                
                
                %obj.get
                
                
                % update time
                obj.grid(:,:,3) = obj.grid(:,:,3) + dt;
                
                
                here = 1;
                obj.stateT('fold');
                
                
                
%                figure(h2);
%                plot(LV)
%                figure(h3);
%                plot(what)
                
                drawnow
                
            end
            here = 1;
        end
        
        
        
        function [] = stateT(obj,finalState)
            
            [finalState] = obj.handleStateTransitionVector(finalState);
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if want to end at folded state
            if strcmp(finalState{1},'fold')
                if strcmp(obj.state{1},'fold')
                    % do nothing
                else
                    % transistion from unfolded - > folded
                    obj.fold();
                end
            else
                if strcmp(obj.state{1},'unfold')
                    % do nothing
                else
                    % transistion from folded -> unfolded
                     obj.unfold();
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.state{1} = finalState{1};
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if want to end at folded state
            if strcmp(finalState{2},'base')
                if strcmp(obj.state{2},'base')
                    % do nothing
                else
                    if strcmp(obj.state{1},'fold');obj.unfold();end
                    % transistion from apex - > base
                    obj.base();
                    if strcmp(obj.state{1},'fold');obj.fold();end
                end
            else
                if strcmp(obj.state{2},'apex')
                    % do nothing
                else
                    % transistion from apex -> base
                    if strcmp(obj.state{1},'fold');obj.unfold();end
                    % transistion from base - > apex
                    obj.apex();
                    if strcmp(obj.state{1},'fold');obj.fold();end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.state{2} = finalState{2};
            
        end
        
        
        function [v] = getVelocity(obj,frame)
            if nargin == 1
                frame = 'apex';
            end
            initState = obj.state;
            obj.stateT({'apex','fold'});
            v = obj.vFunc(obj.grid);
            obj.stateT(initState);
            v = reshape(v,obj.sz(1:2));
            
            % this needs to be fixed to max later -- FIX
            if strcmp(frame,'base')
                v = bsxfun(@plus,-v,v(:,1));
            end
        end
        
        
        function [] = render(obj,h,toQuiver)
            initState = obj.state;
            obj.stateT('unfold');

            x = obj.grid(:,:,5);
            y = obj.grid(:,:,6);
           
            %{
            TAN_L = (dxdl.^2 + dydl.^2).^-.5;
            dxdl = dxdl.*TAN_L;
            dydl = dydl.*TAN_L;


            TAN_W = (dxdw.^2 + dydw.^2).^-.5;
            dxdw = dxdw.*TAN_W;
            dydw = dydw.*TAN_W;

            %}
            figure(h);
            skip_w = 10;
            for e = 1:skip_w:size(obj.grid,1)
                plot(obj.grid(e,:,5),obj.grid(e,:,6),'c')
                hold on
            end
            plot(obj.grid(1,:,5),obj.grid(1,:,6),'b')

            skip_l = 10;
            for e = 1:skip_l:size(obj.grid,2)
                plot(obj.grid(:,e,5),obj.grid(:,e,6),'m')
            end
            plot(obj.grid(:,1,5),obj.grid(:,1,6),'r')
            
            if toQuiver
                [F] = getFrame(obj);
                quiver(obj.grid(1:skip_w:end,1:skip_l:end,5),...
                       obj.grid(1:skip_w:end,1:skip_l:end,6),...
                       F(1:skip_w:end,1:skip_l:end,1,1),...
                       F(1:skip_w:end,1:skip_l:end,2,1),'k')
                %quiver(obj.grid(:,:,5),obj.grid(:,:,6),dydw,dxdw,'y')
            end
            
            %{
            quiver(obj.grid(:,:,1),obj.grid(:,:,2),dydl,dxdl,'k')


            quiver(obj.grid(:,:,1),obj.grid(:,:,2),dydw,dxdw,'y')
            %}


            %{
            plot(vT(2,3),vT(1,3),'g.')
            quiver(vT(2,3),vT(1,3),vT(2,1),vT(1,1),10,'Color','g')

            %plot(0,0,'g*');

            %axis([-5 5 0 100]);
            hold off
            
            %}
            axis equal
            drawnow
            %waitforbuttonpress
            obj.stateT(initState);
            hold off
        end
        
        
        
        function [measurement] = getMeasurement(obj,type,measurement,h)
            if nargin <= 2
                measurement.cell_profile = [];
                measurement.midline_length = [];
                measurement.tip_angle = [];
                measurement.midline_position = [];
                measurement.velocity = [];
            end
            
            if nargin >= 2
                if isempty(measurement)
                    measurement.cell_profile = [];
                    measurement.midline_length = [];
                    measurement.tip_angle = [];
                    measurement.midline_position = [];
                    measurement.velocity = [];
                end
            end
            
            if nargin == 3
                h = [];
            end
            
            initS = obj.state;
            obj.stateT('unfold');
            
            switch type
                case 'cell_profile'
                    cellStrip = obj.grid((end-1)/2,:,5:6);
                    dx = gradient(cellStrip(:,:,1));
                    dy = gradient(cellStrip(:,:,2));
                    m = (dx.^2 + dy.^2).^.5;
                    measurement.cell_profile = [measurement.cell_profile;m];
                case 'midline_length'
                    m = obj.grid((end-1)/2,end,1);
                    measurement.midline_length = [measurement.midline_length;m];
                case 'tip_angle'
                    m = obj.grid((end-1)/2,end,4);
                    measurement.tip_angle = [measurement.tip_angle;m];
                case 'width_profile'
                    %measurement = gradient(
                case 'midline_position'
                    m = obj.grid((end-1)/2,:,[5 6]);
                    m = permute(m,[3 2 1]);
                    measurement.midline_position = cat(3,measurement.midline_position,m);
                case 'velocity'
                    midline = permute(measurement.midline_position,[1 3 2]);
                    tc = squeeze(bsxfun(@minus,midline,midline(:,:,end)));
                    dis = squeeze(sum(tc.*tc,1).^.5);
                    vel = diff(dis,1,1);
                    dis(end,:) = [];
                    measurement.velocity = cat(3,dis,vel);
                case 'display'
                    measurement = obj.getMeasurement('cell_profile',measurement,h);
                    measurement = obj.getMeasurement('midline_length',measurement,h);
                    measurement = obj.getMeasurement('tip_angle',measurement,h);
                    measurement = obj.getMeasurement('midline_position',measurement,h);
                    measurement = obj.getMeasurement('velocity',measurement,h);
                    
                    
                    
                    % perspective of the length to base
                    av = obj.getVelocity('apex');
                    initState = obj.state;
                    obj.stateT('apex');
                    mid = obj.grid((end-1)/2,:,1);
                    obj.stateT(initState);
                    
                    
                    
                    figure(h);
                    subplot(2,2,1)
                    plot(measurement.cell_profile(end,:));
                    subplot(2,2,2);
                    plot(measurement.midline_length)
                    subplot(2,2,3)
                    plot(measurement.tip_angle)
                    subplot(2,2,4)
                    dis = measurement.velocity(:,:,1);
                    vel = measurement.velocity(:,:,2);
                    plot(dis(:),vel(:)*100,'.')
                    hold on
                    %plot(mid,av,'r')
                    hold off
                    drawnow
                    
            end
            obj.stateT(initS);
                   
        end
        
        function [F,DX] = getFrame(obj)
            initState = obj.state;
            obj.stateT('unfold');

            x = obj.grid(:,:,5);
            y = obj.grid(:,:,6);
            [dxdl,dxdw] = gradient(x);
            [dydl,dydw] = gradient(y);

            dL = (dxdl.^2 + dydl.^2).^.5;
            
            dxdl = dxdl.*dL.^-1;
            dydl = dydl.*dL.^-1;


            dW = (dxdw.^2 + dydw.^2).^.5;
            dxdw = dxdw.*dW.^-1;
            dydw = dydw.*dW.^-1;
            obj.stateT(initState);
            
            
            DX = cat(3,dL,dW);
            F = cat(4,cat(3,dxdl,dydl),cat(3,dxdw,dydw));
        end
    end
    
    
    
    methods (Access = private)
        
        
        function [L] = calcLengthFromGrid(obj)
            initState = obj.state;
            obj.stateT('unfold');

            x = obj.grid(:,:,5);
            y = obj.grid(:,:,6);
            [dxdl,dxdw] = gradient(x);
            [dydl,dydw] = gradient(y);
            
            L = cumsum((dxdl.^2 + dydl.^2).^.5,2);
        end
        
        function [] = fold(obj)
            obj.grid = reshape(obj.grid,[prod(obj.sz(1:2)) obj.sz(end)]);
        end
        
        function [] = unfold(obj)
            obj.grid = reshape(obj.grid,obj.sz);
        end
        
        function [] = base(obj)
            obj.grid(:,:,1) = bsxfun(@plus,-obj.grid(:,:,1),obj.grid(:,1,1));
        end
        
        function [] = apex(obj)
            obj.grid(:,:,1) = bsxfun(@plus,-obj.grid(:,:,1),obj.grid(:,end,1));
        end
        
        function [stateVector] = handleStateTransitionVector(obj,stateVector)
            
            if ~iscell(stateVector)
                tmp = stateVector;
                stateVector = {};
                stateVector{1} = tmp;
            end
            
            for e = 1:numel(stateVector)
                if strcmp(stateVector{e},'fold') || strcmp(stateVector{e},'unfold')
                    TransitionType{e} = 'foldType';
                elseif strcmp(stateVector{e},'base') || strcmp(stateVector{e},'apex')
                    TransitionType{e} = 'centerType';
                else
                    fprintf(['Transition Error']);
                end
            end
            
            
            if numel(TransitionType) == 1
                if strcmp(TransitionType{1},'foldType')
                    TransitionType{2} = 'centerType';
                    stateVector{2} = obj.state{2};
                else
                    TransitionType{2} = 'foldType';
                    stateVector{2} = obj.state{1};
                end
            end
            
            
            if strcmp(TransitionType{1},'centerType')
                stateVector = flip(stateVector,2);
            end
            
            
            
            
        end
    end
    
end