
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dig for the raw data for extracting the wave front data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['*****************************************************.\n'])
fprintf(['Looking for wave front data.\n'])
matFilePath = '/mnt/snapper/nate/AminoAcidWave/Extraction/';
matFileList = {};
FileExt = {'mat'};
matFileList = sdig(matFilePath,matFileList,FileExt,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['*****************************************************.\n'])
fprintf(['Looking for wave front data.\n'])
fprintf(['*****************************************************.\n'])
fprintf(['Loading wave front data.\n']);
randN = 50;
D = [];
for s = 1:numel(matFileList)
    tmpList = matFileList{s};
    tmpList = tmpList(randperm(numel(tmpList)));
    for e = 1:randN
        try
            a = load(tmpList{e});
            D = cat(4,D,a.rise_data_out);
            fprintf(['Done loading:' num2str(e) ':' num2str(randN) ':' num2str(s) ':' num2str(numel(matFileList)) '\n'])
        catch
        end
    end
end
fprintf(['*****************************************************.\n'])
fprintf(['Loading wave front data.\n']);
fprintf(['*****************************************************.\n'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make and remove outliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tr = 1;
reduce1 = 3;
reduce2 = 3;
tmpQ1 = D;
fprintf(['*****************************************************.\n'])
fprintf(['Building model and removing outlier.\n']);
fprintf(['*****************************************************.\n'])
while tr ~= size(tmpQ1,4)
    tr = size(tmpQ1,4);
    qSZ = size(tmpQ1);
    tmpQ1 = reshape(tmpQ1,[prod(qSZ(1:2)) prod(qSZ(3:4))]);
    [U_rise1,E_rise1,L_rise1] = PCA_FIT_FULL_Tws(tmpQ1,reduce1);
    C_rise1 = PCA_REPROJ_T(tmpQ1,E_rise1,U_rise1);
    C_rise1 = reshape(C_rise1,[reduce1 qSZ(3:4)]);
    qSZ1 = size(C_rise1);
    tmpC1 = reshape(C_rise1,[prod(qSZ1(1:2)) prod(qSZ1(3))]);
    [U_rise2,E_rise2,L_rise2] = PCA_FIT_FULL_Tws(tmpC1,reduce2);
    C_rise2 = PCA_REPROJ_T(tmpC1,E_rise2,U_rise2);
    C_rise2 = reshape(C_rise2,[reduce2 qSZ1(3)]);
    TF = isoutlier(C_rise2(1,:));
    tmpQ1 = reshape(tmpQ1,qSZ);
    tmpQ1(:,:,:,TF) = [];
end
fprintf(['*****************************************************.\n'])
fprintf(['Building model and removing outlier.\n']);
fprintf(['*****************************************************.\n'])

%% look at level
close all

ep = .01;
for tr =1:size(tmpQ1,4)
    for t = 1:size(tmpQ1,3)
        value = ba_interp2(tmpQ1(:,:,t,tr),25.5,25.5);
        msk = tmpQ1(:,:,t,tr) > (value - ep) & tmpQ1(:,:,t,tr) < (value + ep);
        imshow(msk,[]);
        drawnow
    end
end
%%

%% mask only - try one pass
tmpQ1 = nQ1;
close all
qSZ = size(tmpQ1);
reduce1 = 3;
tmpQ1 = reshape(tmpQ1,[prod(qSZ(1:2)) prod(qSZ(3:4))]);
[U_rise1,E_rise1,L_rise1] = PCA_FIT_FULL_Tws(tmpQ1,reduce1);
C_rise1 = PCA_REPROJ_T(tmpQ1,E_rise1,U_rise1);
C_rise1 = reshape(C_rise1,[reduce1 qSZ(3:4)]);
qSZ1 = size(C_rise1);
reduce2 = 3;
tmpC1 = reshape(C_rise1,[prod(qSZ1(1:2)) prod(qSZ1(3))]);
[U_rise2,E_rise2,L_rise2] = PCA_FIT_FULL_Tws(tmpC1,reduce2);
C_rise2 = PCA_REPROJ_T(tmpC1,E_rise2,U_rise2);
C_rise2 = reshape(C_rise2,[reduce2 qSZ1(3)]);
TF = isoutlier(C_rise2(1,:));
tmpQ1 = reshape(tmpQ1,qSZ);
tmpQ1(:,:,:,TF) = [];
plot3(C_rise2(1,:),C_rise2(2,:),C_rise2(3,:),'.')
