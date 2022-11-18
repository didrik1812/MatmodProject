classdef DiffusionComponentInputParamsGlia < ComponentInputParams

    properties
        
        D
        
    end
    
    
    methods

        function paramobj = DiffusionComponentInputParamsGlia(jsonstruct)
            paramobj = paramobj@ComponentInputParams(jsonstruct);
        end
        
    end
    

    
end
