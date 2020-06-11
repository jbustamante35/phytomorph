classdef htMap
    
    
    
    
    methods (Static)
        
        function [] = getXferCommand(fileName)
            location = getLocation(FilePath);
            fprintf(['search@' location '@' FilePath '\n']);
            switch location
                case 'irods'
                   
                case 'cold'
                    
                case 'local'
                    
                case 'chtc'
                    
                case 'squid'
                    
            end
        end
        
        
    end
end