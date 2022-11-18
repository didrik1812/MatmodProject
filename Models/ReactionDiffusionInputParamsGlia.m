classdef ReactionDiffusionInputParamsGlia < InputParams

    properties
        
        G
        
        T % parameter for diffusion component model T
        N % parameter for diffusion component model N
        T_N % parameter for diffusion component model T_N
        N_I % parameter for diffusion component model N_I
        
        k_1 % reaction rate k_1
        k_2 % reation rate k_{-1}
        k_3 % reaction rate k_2

    end
    
    methods
        
        function paramobj = ReactionDiffusionInputParamsGlia(jsonstruct)
            
            paramobj = paramobj@InputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);
            
            paramobj.T = DiffusionComponentInputParamsGlia(pick('T'));
            paramobj.N = DiffusionComponentInputParamsGlia(pick('N'));
            paramobj.T_N = DiffusionComponentInputParamsGlia(pick('T_N'));
            paramobj.N_I = DiffusionComponentInputParamsGlia(pick('N_I'));

            paramobj = paramobj.validateInputParams();
            
        end

        function paramobj = validateInputParams(paramobj)

            if ~isempty(paramobj.G)
                paramobj.T.G = paramobj.G;
                paramobj.N.G = paramobj.G;
                paramobj.T_N.G = paramobj.G;
                paramobj.N_I.G = paramobj.G;
            end
            
        end
        
    end
    
end