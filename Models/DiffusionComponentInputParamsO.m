classdef DiffusionComponentInputParamsO < ComponentInputParams

    properties
        
        D
        
    end
    
    
    methods

        function paramobj = DiffusionComponentInputParamsO(jsonstruct)
            paramobj = paramobj@ComponentInputParams(jsonstruct);
        end
        
    end
    

    
end
