function [] = mainTRiP(FileList,disp,oPath)
    myAutoCropper_TRip(FileList,disp,oPath);
    estimateAll(oPath,oPath);
    modelFitAll();
end

%{
    FilePath = 'W:\';
    FilePath = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/TRiP/TRiP_1.0/input/';
    FileList = {};
    FileExt = {'tif','TIF','jpg'};
    FileList = gdig(FilePath,FileList,FileExt,1);
    mainTRiP(FileList,false,'./moutput/');
%}


