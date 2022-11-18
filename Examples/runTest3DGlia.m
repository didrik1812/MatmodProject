close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('DiffusionGlia.json');
jsonstruct = jsondecode(jsonfile);

%Set N.D to around 70 in jsonfile for 2D and 100 to 3D

paramobj       = ReactionDiffusionInputParamsGlia(jsonstruct);
paramobj.k_1   = paramobj.k_1;
paramobj.k_2   = paramobj.k_2;
paramobj.k_3   = paramobj.k_3;
paramobj.T.D   = paramobj.T.D;
paramobj.N.D   = paramobj.N.D;
paramobj.T_N.D   = paramobj.T_N.D;
paramobj.N_I.D = paramobj.N_I.D;

G = Cylindergrid();
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusionGlia(paramobj);

% setup schedule
total = 3e-5;
n  = 10;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);


% setup initial state
nc     = G.cells.num;
%disp(nc)
vols   = G.cells.volumes;
initCT = 1000*((meter)^2)/sum(G.cells.volumes(1:1218)); 
initCT_2D=1000/sum(G.cells.volumes(1:3:1218));

V_layer_bottom = sum(G.cells.volumes(756:757));
V_layer_bottom_middle_cell=G.cells.volumes(756:757);
V_2D=G.cells.volumes(829:830);

m=1;
initCN=5000*m;%fix this according to dimension
initcase = 1;

switch initcase
  case 1
    cT            = zeros(nc, 1);
    cT(1:1218)     = initCR;
    cN            = zeros(nc, 1);
    cN(10963:12180) = initCN/V_layer_bottom;
    cT_N          = zeros(nc, 1);
    cN_I          = zeros(nc, 1);
  case 2
    cT            = zeros(nc, 1);
    cT(1:1218)     = initCT;
    cN            = zeros(nc, 1);
    cN(10968:10969) = initCN/V_layer_bottom_middle_cell;
    cN_I          = zeros(nc, 1);
  case 3
    cT            = zeros(nc, 1);
    cT(1:3:1218)     = initCT;
    cN            = zeros(nc, 1);
    cN(829:830) = initCN/V_2D;
    cN_I          = zeros(nc, 1);
        
end


initstate.T.c   = cT;
initstate.N.c   = cN;
initstate.T_N.c = cT_N;
initstate.N_I.c = cN_I;

G = model.G;
receptorCells = (1 : 100);
injectionCells = (66:66);
% plot injection and receptor cells
%[16,20,66,74] %Some indexes of cells on the border
figure
%plotToolbar(G,states)
plotGrid(G, 'facecolor', 'none');
plotGrid(G, injectionCells, 'facecolor', 'red');
%plotGrid(G, receptorCells, 'facecolor', 'yellow');
view(33, 26);
return
% run simulation

nls = NonLinearSolver();
nls.errorOnFailure = false;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls);



%%

% Remove empty states (could have been created if solver did not converge)
ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

figure(1); figure(2); figure(3); figure(4);

for istate = 1 : numel(states)

    state = states{istate};

    set(0, 'currentfigure', 1);
    cla
    plotCellData(model.G, state.T.c);view(30,60);
    colorbar
    title('T concentration')
    
    set(0, 'currentfigure', 2);
    cla
    plotCellData(model.G, state.N.c); view(30,60);
    colorbar
    title('N concentration')

    set(0, 'currentfigure', 3);
    cla
    plotCellData(model.G, state.N_I.c);view(30,60);
    colorbar
    title('N_I concentration')
    
    set(0, 'currentfigure', 4);
    cla
    plotCellData(model.G, state.T_N.c);view(30,60);
    colorbar
    title('T_N concentration')

    drawnow
    pause(0.1);
    
end
