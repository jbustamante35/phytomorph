%%
FilePath = '/mnt/tetra/nate/MATforshiheng/';
FileList = {};
FileExt = {'mat'};
tic
FileList = gdig(FilePath,FileList,FileExt,1);
toc
%% this was downloaded and run 
FilePath = '/mnt/tetra/nate/MATforshiheng/';
FileList = {};
FileExt = {'mat'};
tic
FileList = gdig(FilePath,FileList,FileExt,1);
toc
%%
for e = 1:numel(FileList)
    %FileList{e}
    if contains(FileList{e},'Bd21-3')
        
        e
    end
end
%%
a = load(FileList{317});
%%
a.out
%%
for ee = 1:numel(FileList)
    
    a = load(FileList{ee});
    [pth,nm,ext] = fileparts(FileList{ee});
    
    KUR = [];
    close all
    h1 = figure;
    h2 = figure;
    SNIP = 21;
    rootDisp = false;
    for tm = 1:numel(a.out)
        try

            for root = 1:numel(a.out{tm}.midlines)
                if rootDisp
                    figure(h1);
                    plot(a.out{tm}.midlines(root).data(1,:),a.out{tm}.midlines(root).data(2,:),'r')
                    hold on
                end


                para = {5};
                K = cwtK_root(a.out{tm}.midlines(root).data',para);

                if tm == 1
                    KUR(:,tm,root) = K.K;
                else
                    KUR(:,tm,root) = K.K(1:size(KUR,1));
                    %{
                    figure(h2);
                    mesh(KUR(SNIP:end-SNIP,:,1));
                    view([0 90]);
                    %}
                end

            end
            if rootDisp
            figure(h1);
            for root = 1:numel(a.out{tm}.midlines)
                plot(a.out{tm}.contours(root).data(1,:),a.out{tm}.contours(root).data(2,:),'b')
                hold on
            end
            end

           %{
            hold off
            drawnow
            %}
        catch
        end

        tm
    end
    try
    %%
    oPath = '/mnt/tetra/nate/curvatureReturn/';
    close all
    SELECT = 300;
    for root = 1:size(KUR,3)
        ksig = -KUR(SNIP:end-SNIP,:,root);
        ksig = imfilter(ksig,fspecial('disk',[5]),'replicate');
        figure;
        mesh(ksig(1:SELECT,:));
        view([0 90]);
        axis([0 size(ksig,2) 0 300 -.006 .006])
        caxis([-.006 .006])
        view([0 90]);
        colorbar
        saveas(gca,[oPath nm '__' num2str(root) '.tif']);
        close all
    end
    catch
    end
end
                
               