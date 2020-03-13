classdef cotModel < matlab.mixin.Heterogeneous

    properties
        fileName;

        tVec;

        alignedImage;
        alignedMask;
        labelImage;

        holesImage;
        numberHoles;

        veinNetwork;

    end

    methods

        function [obj] = cotModel(fileName,tVec)
            obj.fileName = fileName;
            obj.tVec = tVec;
        end

        function [] = decorateVeins(obj)
            obj.veinNetwork.decorateGraph();
        end

        function [] = view(obj)
            renderSingleCotImage(obj,[])
            
            %imshow(obj.alignedImage,[]);
            %drawnow

        end
        

    end

end