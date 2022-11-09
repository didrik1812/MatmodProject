close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('diffusion2.json');
jsonstruct = jsondecode(jsonfile);

paramobj = ReactionDiffusionInputParams(jsonstruct);

G = Cylindergrid();
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusion(paramobj);


% setup schedule
total = 100;
n  = 100;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

% setup initial state

nc = G.cells.num;
vols = G.cells.volumes;

initcase = 1;
switch initcase
  case 1
    cR      = zeros(nc, 1);
    cR(1:578)   = sum(vols);
    cN      = zeros(nc, 1);
    cN(5202:5780) = sum(vols);
    cR_N = zeros(nc, 1);
  case 2
    cR = ones(nc, 1);
    cN = ones(nc, 1);
    cR_N = zeros(nc, 1);
end

initstate.R.c = cR;
initstate.N.c = cN;
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
