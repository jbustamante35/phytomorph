classdef compressionLayer < configurableLayer
    
    properties
        referenceLayer;
        referenceGrade;
        compressLevel = 5;
        U;
        E;
        L;
    end
    
    methods


        function [obj] = compressionLayer(layerToCompress)
            name = ['compress-' layerToCompress.name];
            obj = obj@configurableLayer(name);
            obj.referenceLayer = layerToCompress;

        end
        
        function [out] = compute(obj,varargin)

            if obj.isConfigured
                tmp = thawTensor(varargin{1},obj.referenceGrade);
                tmp = freezeTensor(tmp);
                res = PCA_REPROJ_T(tmp,obj.E,obj.U,false);
            end
            out = varargin;
            out{1} = res;
        end

        function [] = setReferenceGrade(obj,rGrade)
            obj.referenceGrade = rGrade;
        end
        
        function [] = setCompressLevel(obj,compressLevel)
            obj.compressLevel = compressLevel;
        end


        function [] = configure(obj,varargin)
            if ~obj.isConfigured

                %{
                flowData = traceNat('bug1');


                flowData.flowDirection = 'r';
                flowData.stopData.targetQueue = obj.topCE;


                flowData.persistData.toPersist = false;
                flowData.persistData.userName = 'nmiller';
                flowData.persistData.location = 'irods';

                flowData.jobData.uid = '1';
                %}


                %{
                flowData.flowDirection = 'r';
                flowData.history = [];

                flowData.stopData.uid = obj.topCE.uid;
                flowData.stopData.goFlag = false;

                flowData.persistData.toPersist = false;
                flowData.persistData.userName = 'nmiller';
                flowData.persistData.location = 'irods';

                flowData.jobData.uid = '1';
                %}

               

                for e = 1:numel(varargin{1})
                    %flowData.flowDirection = 'r';


                    tmpFluid = digitalFluid('tmpFluid',varargin{1}{e});


                    tmp = tmpFluid.flow(obj.topCE,'r');


                    %{
                    tmpFluid.flowDirection = 'r';
                    tmpFluid.stopData.targetQueue = obj.topCE;
                    tmp = obj.flow(tmpFluid);
                    %}
                    tmp = tmp.getCurrent();
                    tmp = tmp{1};



                    tmp = thawTensor(tmp,obj.referenceGrade);
                    tmp = freezeTensor(tmp);
                    if e == 1
                        dataStack = zeros(size(tmp,1),numel(varargin{1}));
                    end
                    dataStack(:,e) = tmp;
                end

                [obj.U,obj.E,obj.L] = PCA_FIT_FULL_Tws(dataStack,obj.compressLevel);
                obj.isConfigured = true;
            end
        end
        
    end
end