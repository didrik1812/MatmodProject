close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('diffusion2.json');
jsonstruct = jsondecode(jsonfile);

%Set N.D to around 70 in jsonfile for 2D and 100 to 3D

paramobj       = ReactionDiffusionInputParams([jsonstruct]);
paramobj.k_1   = paramobj.k_1;
paramobj.k_2   = paramobj.k_2;
paramobj.N.D   = paramobj.N.D;
paramobj.R.D   = paramobj.R.D;
paramobj.R_N.D = paramobj.R_N.D;

G = Cylindergrid();
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusion(paramobj);

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
initCR = 1000*((meter)^2)/sum(G.cells.volumes(1:1218)); 
initCR_2D=1000/sum(G.cells.volumes(1:3:1218));
V_layer_bottom = sum(G.cells.volumes(756:757));
V_layer_bottom_middle_cell=G.cells.volumes(756:757);
V_2D=G.cells.volumes(829:830);
m=1;
initCN=5000*m;%fix this according to dimension
initcase = 3;

switch initcase
  case 1
    cR            = zeros(nc, 1);
    cR(1:1218)     = initCR;
    cN            = zeros(nc, 1);
    cN(10963:12180) = initCN/V_layer_bottom;
    cR_N          = zeros(nc, 1);
  case 2
    cR            = zeros(nc, 1);
    cR(1:1218)     = initCR;
    cN            = zeros(nc, 1);
    cN(10968:10969) = initCN/V_layer_bottom_middle_cell;
    cR_N          = zeros(nc, 1);
  case 3
    cR            = zeros(nc, 1);
    cR(1:3:1218)     = initCR;
    cN            = zeros(nc, 1);
    cN(829:830) = initCN/V_2D;
    cR_N          = zeros(nc, 1);
        
end


initstate.R.c   = cR;
initstate.N.c   = cN;
initstate.R_N.c = cR_N;

%G = model.G;
%receptorCells = (1 : 100);
%injectionCells = (829:830);
% plot injection and receptor cells

%figure
%plotToolbar(G,states)
%plotGrid(G, 'facecolor', 'none');
%plotGrid(G, injectionCells, 'facecolor', 'red');
%plotGrid(G, receptorCells, 'facecolor', 'yellow');
%view(33, 26);
%return
% run simulation

nls = NonLinearSolver();
nls.errorOnFailure = false;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls);



%%

% Remove empty states (could have been created if solver did not converge)
ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

figure(1); figure(2); figure(3);

for istate = 1 : numel(states)

    state = states{istate};

    set(0, 'currentfigure', 1);
    cla
    plotCellData(model.G, state.R.c);view(30,60);
    colorbar
    title('R concentration')
    
    set(0, 'currentfigure', 2);
    cla
    plotCellData(model.G, state.N.c); view(30,60);
    colorbar
    title('N concentration')

    set(0, 'currentfigure', 3);
    cla
    plotCellData(model.G, state.R_N.c);view(30,60);
    colorbar
    title('R-N concentration')

    drawnow
    pause(0.1);
    
end
