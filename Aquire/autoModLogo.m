
I = double(imread('~/phytoMorphTK/pmLogo.png'))/255;
J = double(imread('~/phytoMorphTK/htCondor.png'))/255;

sz = size(J);
J = reshape(J,[prod(sz(1:2)) sz(3)]);
fidx = find(all(bsxfun(@eq,J,[1 1 1]),2));
J(fidx,:) = repmat([.94 .94 .94],[numel(fidx) 1]);
J = reshape(J,sz);

sz = size(I);
I = reshape(I,[prod(sz(1:2)) sz(3)]);
fidx = find(all(bsxfun(@eq,I,[1 1 1]),2));
I(fidx,:) = repmat([.94 .94 .94],[numel(fidx) 1]);
I = reshape(I,sz);

%%
close all
h = size(J,1);
w = h/size(I,1)*size(I,2);
b = round((size(J,2) - w)/2);
buffer = .94*ones(h,b,3);
Imod = imresize(I,[h w]);
top = [buffer,Imod,buffer];
modI = [J;top];
imshow(modI,[]);
drawnow
imwrite(modI,['~/phytoMorphTK/logo.png']);
%%
close all
[K,map,alpha] = imread('~/phytoMorphTK/cyverseLogo.png');
map(1,:) = [.94 .94 .94];
K = ind2rgb(K,map);
imwrite(K,'~/phytoMorphTK/cyverse.png');