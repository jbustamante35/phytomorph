function [xy] = findHorizontalLines(I)


    gr = rgb2gray(I);
    gr = imfilter(gr,fspecial('disk',11),'replicate');
    BW = edge(gr);

    
    
    
    [H,T,R] = hough(BW','Theta',linspace(-10,10,100));
    P  = houghpeaks(H,1);
    lines = houghlines(BW',T,R,P,'FillGap',400,'MinLength',100);
    %{
    for k = 1:numel(lines)
        xy = [lines(k).point1; lines(k).point2];
        plot(xy(:,2),xy(:,1),'LineWidth',2,'Color','green');
        
    end
    %}
    
    xy1 = [lines(1).point1; lines(1).point2];
    CHOP = mean(xy1(:,1));
    mBW = BW(1:(CHOP-10),:);
    
    [H,T,R] = hough(mBW','Theta',lines.theta);
    P  = houghpeaks(H,1);
    lines = houghlines(mBW',T,R,P,'FillGap',400,'MinLength',10);
    %xy2 = [lines(1).point1; lines(1).point2];
    xy2 = [];
    for k = 1:numel(lines)
        xy2 = [lines(k).point1; lines(k).point2];
        %plot(xy(:,2),xy(:,1),'LineWidth',2,'Color','red');
        
    end
    
    xy = cat(3,xy1,xy2);
    
end