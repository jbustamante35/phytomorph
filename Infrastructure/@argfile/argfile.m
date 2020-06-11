classdef argfile < file
    properties
        numArgs;
    end
    
    methods
        function [this] = argfile(fName,numArgs)
            this@file(fName);
            this.numArgs = numArgs;
        end
        function [args] = load(this)
            a = load([this.filePath this.fileName],'in');
            args = a.in;
        end
    end
end