function [I] = makeLabelText(text,sz,fontSize)
    cmd = ['convert -size ' num2str(sz(1)) 'x' num2str(sz(2)) ' -background white -pointsize ' num2str(fontSize) ' -fill black -gravity West caption:''' ...
        text ''' ~/test.tif'];
    system(cmd);
    I = imread('~/test.tif');
    I = double(I);
    I = bindVec(I(:,:,1));
end