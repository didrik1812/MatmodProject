close all

mrstModule add ad-core mrst-gui 

%jsonfile = fileread('DiffusionGlia.json');
%jsonstruct = jsondecode(jsonfile);

paramobj       = ReactionDiffusionInputParamsGlia([]);
paramobj.k_1   = 4e6*(mol/litre)*(1/second);
paramobj.k_2   = 5*(1/second);
paramobj.k_3   = 5e10*(mol/litre)*(1/second); %Find units here
paramobj.N.D   = 8e-7*(meter^2/second);
paramobj.T.D   = 0*(meter^2/second);
paramobj.T_N.D = 0*(meter^2/second);
paramobj.N_I.D = 0*(meter^2/second);

G = Cylindergrid();
G = computeGeometry(G);

max_volume = max(G.cells.volumes);

edgeCells = [];
for i = 1:G.cells.num
    if G.cells.volumes(i) >= max_volume - max_volume * 1/100
        edgeCells(end + 1) = i;
    end
end

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusionGlia(paramobj);

% setup schedule
total = 1*nano*second;
n     = 100;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

A = 6.02214076e23; % Avogadro constant

gliaCells = (edgeCells);
injectionCells = (1 : 12180); 

% setup initial state
nc     = G.cells.num;
%disp(nc)
vols   = G.cells.volumes;
initCT = (1000/A)*((meter)^2)/sum(G.cells.volumes(gliaCells)); 

V      = sum(G.cells.volumes(injectionCells));
initCN = (5000/A)/V;
initcase = 1;

switch initcase
  case 1
    cT            = zeros(nc, 1);
    cT(gliaCells)     = initCT;
    cN            = zeros(nc, 1);
    cN(injectionCells) = initCN;
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

nls = NonLinearSolver();
nls.errorOnFailure = false;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls);

figure
plotToolbar(G,states)
return

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
