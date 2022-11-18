classdef ReactionDiffusionGlia < BaseModel

    properties
        
        T % Diffusion model for component R
        N % Diffusion model for component N
        T_N %Diffuison model for component T_N
        N_I % Diffusion model for result N_I
        
        k_1 % reaction constant k_1
        k_2 % reaction constant k_2
        k_3 % reaction constant k_3
        
    end
    
    methods
        
        function model = ReactionDiffusionGlia(paramobj)
            
            model = model@BaseModel();
            
            % All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            %% Setup the model using the input parameters
            fdnames = {'G'            , ...
                       'k_1'          , ...
                       'k_2'          , ...
                       'k_3'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.T = DiffusionComponentGlia(paramobj.T);
            model.N = DiffusionComponentGlia(paramobj.N);
            model.T_N = DiffusionComponentGlia(paramobj.T_N);
            model.N_I = DiffusionComponentGlia(paramobj.N_I);
            
        end

        
        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            %% Temperature dispatch functions
            fn = @ReactionDiffusionGlia.updateSourceTerm;
            
            inputnames = {{'T', 'c'}, ...
                          {'N', 'c'}, ...
                          {'T_N', 'c'}, ...
                          {'N_I', 'c'}};
            model = model.registerPropFunction({{'T', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'N', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'T_N', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'N_I', 'source'} , fn, inputnames});
            
        end

        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@PhysicalModel(model);
            forces.none = [];
        end
        
        function state = updateSourceTerm(model, state)

            k_1 = model.k_1;
            k_2 = model.k_2;
            k_3 = model.k_3;
            vols = model.G.cells.volumes;
            
            cT = state.R.c;
            cN = state.N.c;
            cT_N = state.T_N.c;
            cN_I = state.N_I.c;

            Re1 = k_1.*vols.*cT.*cN;
            Re2 = k_2.*vols.*cT_N;
            Re3 = k_3.*vols.*cT_N;

            state.T.source = -Re1 +Re2;
            state.N.source = -Re1 +Re2;
            state.T_N.source = Re1 -Re2 - Re3;
            state.N_I.source = Re3;
            
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            state = model.initStateAD(state);

            state.T = model.T.updateFlux(state.T);
            state.N = model.N.updateFlux(state.N);
            state.T_N = model.T_N.updateFlux(state.T_N);
            state.N_I = model.N_I.updateFlux(state.N_I);
            
            state = model.updateSourceTerm(state);

            state.T = model.T.updateMassAccum(state.T, state0.T, dt);
            state.N = model.N.updateMassAccum(state.N, state0.N, dt);
            state.T_N = model.T_N.updateMassAccum(state.T_N, state0.T_N, dt);
            state.N_I = model.N_I.updateMassAccum(state.N_I, state0.N_I, dt);
            
            state.T = model.T.updateMassConservation(state.T);
            state.N = model.N.updateMassConservation(state.N);
            state.T_N = model.T_N.updateMassConservation(state.T_N);
            state.N_I = model.N_I.updateMassConservation(state.N_I);
            
            eqs = {}; types = {}; names = {};
            
            eqs{end + 1}   = state.T.massCons;
            names{end + 1} = 'massCons T';
            types{end + 1} = 'cell';
            
            eqs{end + 1}   = state.N.massCons;
            names{end + 1} = 'massCons N';
            types{end + 1} = 'cell';
            
            eqs{end + 1}   = state.T_N.massCons;
            names{end + 1} = 'massCons T_N';
            types{end + 1} = 'cell';

            eqs{end + 1}   = state.N_I.massCons;
            names{end + 1} = 'massCons N_I';
            types{end + 1} = 'cell';
                        
            primaryVars = model.getPrimaryVariables();

            %% Setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        
        function state = initStateAD(model, state)
        % initialize a new cleaned-up state with AD variables

            % initStateAD in BaseModel erase all fields
            newstate = initStateAD@BaseModel(model, state);
            newstate.time = state.time;
            state = newstate;
            
        end 

        
        function primaryvarnames = getPrimaryVariables(model)

            primaryvarnames = {{'T', 'c'}, ...
                               {'N', 'c'}, ...
                               {'T_N', 'c'}, ...
                               {'N_I', 'c'}};
            
        end
        

        function model = validateModel(model, varargin)
        % nothing special to do
        end


    end

end