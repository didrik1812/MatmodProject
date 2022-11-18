close all

mrstModule add ad-core mrst-gui 

% jsonfile = fileread('diffusion2.json');
% jsonstruct = jsondecode(jsonfile);

paramobj       = ReactionDiffusionInputParams([]);
paramobj.k_1   = 4e6*(mol/litre)*(1/second);
paramobj.k_2   = 5*(1/second);
paramobj.N.D   = 8e-7*(meter^2/second);
paramobj.R.D   = 0*(meter^2/second);
paramobj.R_N.D = 0*(meter^2/second);

G = Cylindergrid();
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusion(paramobj);

% setup schedule
total = 1*micro*second;
n     = 100;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

G = model.G;
receptorCells = (1 : 1218)';
injectionCells = (10963 : 12180)';


doplot = false;
if doplot
    % plot injection and receptor cells
    figure
    plotGrid(G, 'facecolor', 'none');
    plotGrid(G, injectionCells, 'facecolor', 'yellow');
    plotGrid(G, receptorCells, 'facecolor', 'red');
    view(33, 26);
    return
end

% setup initial state

A = 6.02214076e23; % Avogadro constant

nc     = G.cells.num;
vols   = G.cells.volumes;

initCR = (1000/A)*((micro*meter)^2)/sum(G.cells.volumes(receptorCells));

V      = sum(G.cells.volumes(injectionCells));
initCN = (5000/A)/V;

initcase = 1;

switch initcase
  case 1
    cR                 = zeros(nc, 1);
    cR(injectionCells) = initCR;
    cN                 = zeros(nc, 1);
    cN(receptorCells)  = initCN;
    cR_N               = zeros(nc, 1);
  case 2
    cR   = ones(nc, 1);
    cN   = ones(nc, 1);
    cR_N = zeros(nc, 1);
end

initstate.R.c   = cR;
initstate.N.c   = cN;
initstate.R_N.c = cR_N;

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
