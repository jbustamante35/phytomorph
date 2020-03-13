function [D,dD,gD,sD] = generateSdomain(DMSZ,disp)
    % select initial velocity
    vXminRange = [-1 1];
    vXmin = vXminRange(1) + (vXminRange(2)-(vXminRange(1)))*rand(1);
    % select delta velocity
    deltaRange = [-1 1];
    deltaV = deltaRange(1) + (deltaRange(2)-deltaRange(1))*rand(1);
    % set max values
    vxMax = vXmin + deltaV;
    %fprintf(['Vmax:Vmin' num2str([vXmin,vxMax]) '\n'])
    % select acceleration
    acRange = [.001 2];
    ax = acRange(1) + (acRange(2)-acRange(1))*rand(1);
    % seet break point
    breakPointTX = 0;
    VX = @(y)(vxMax*(1+exp(-ax*(y-breakPointTX))).^-1) + vXmin;
    VX = @(y).5*ones(size(y));
    vyMax = 1;
    VY = @(y)vyMax*ones(size(y));

    [y,x] = ndgrid(linspace(-DMSZ,DMSZ,2*DMSZ+1),linspace(-DMSZ,DMSZ,2*DMSZ+1));



    NX = cumsum(VX(y),1)+x;
    NY = y(1,1)+cumsum(VY(y),1);



    NY = NY - NY(DMSZ+1,DMSZ+1);
    NX = NX - NX(DMSZ+1,DMSZ+1);
    
    %{
    NY = x1;
    NX = x2;
    %}

    [NYpNX,NYpNY] = gradient(NY);
    [NXpNX,NXpNY] = gradient(NX);

    D = [NX(:) NY(:) ones(size(NX(:)))];
    sD = [x(:) y(:) ones(size(NX(:)))];

    %dD = cat(3,cat(4,dx1pdx1,dx1pdx2),cat(4,dx2pdx1,dx2pdx2));
    %dD = cat(3,cat(4,NXpNY,NXpNX),cat(4,NYpNY,NYpNX)); % good one
    dD = cat(4,cat(3,NXpNX,NXpNY),cat(3,NYpNX,NYpNY));
    szdD = size(dD);
    dD = squeeze(reshape(dD,[prod(szdD(1:2)) szdD(3:4)]));
    %dD = permute(dD,[1 3 2]);
    gD = cat(3,NX,NY);


    if disp
        for s = 1:size(NY,1)
            plot(NX(s,:),NY(s,:),'r')
            hold all
        end

        for s = 1:size(NY,1)
            plot(NX(:,s),NY(:,s),'b')
            hold all
        end
        drawnow
    end
end