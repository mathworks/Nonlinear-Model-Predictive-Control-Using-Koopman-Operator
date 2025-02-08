%% A MATLAB Tutorial for Control Koopman Operator

% The Koopman operator is a linear infinite-dimensional operator that can
% represent nonlinear dynamics through its action on some observable
% functions. The space created by the observable functions is sometimes
% called "the lifted space". The Koopman operator works in the lifted space
% where it applies a linear regression to calculate the next step model
% prediction. The key idea is that we can still obtain a linear
% representation of the controlled nonlinear system in a higher-dimensional
% space of observables. To help the reader understand the Koopman lifting
% concept, please check the separate script: Koopman_lifting_concept.m
%
% Here we demonstrate the "Control Koopman Operator" or
% "Input-Output Koopman Operator" concept and apply it to control a nonlinear
% system.
%
% Some practical considerations:
% -Choice of observables is crucial - you need enough observables to capture
% both state and control interactions
%
% - The quality of the Koopman approximation depends on the training data coverage
%
% - For control design, you might want to consider specialized observables
% based on your control objectives
%
% Here we show how to design specialized observables for control
% objectives using a simple nonlinear system where we want to track a
% reference signal. We'll use observables that explicitly incorporate
% tracking error and control effort.
%
% This example demonstrates several key concepts (for the meaning of
% observable obsFcn functions, see code down below):
%
% Special observables for control:
% - Tracking error observables
% - Control effort observables
% - Energy-like observables
% - Reference tracking observables
% - Cross-terms for better dynamics capture
%
% The observables are chosen to:
% - Directly represent control objectives (tracking error)
% - Capture important system properties (energy, stability)
% - Enable efficient control computation
%
% MPC design using Koopman model and nlmpcMultistage object availble from Model-Predictive Control Toolbox:
% - Prediction using learned Koopman operator
% - Cost function incorporating tracking and control objectives

% Copyright 2025, The MathWorks Inc.

%% Collect training data
% Using a nonlinear system under test and a feedback control that keeps
% it stable during data collection experiment.

% Define the time span and discretization
dt = 0.01; % sample time
T = 10; % experiment length
t = 0:dt:T;

% Define the state dynamics for the discrete-time nonlinear system (toy model):
f = @(x,u) [x(1) + dt*(x(2) + 0.5*u);
    x(2) + dt*(-sin(x(1)) - 0.2*x(2) + u)];

% Initial conditions and reference.
% Note that we use an extended system, to add integral action.
x0 = [2; 0];
nx = length(x0);
x0_extended = [x0; 0]; % [ position; velocity; integralOfError ]
nx_extended = length(x0_extended);

ref_position = 1.5*sin(0.5*t');  % Reference trajectory to track

% Simple PD controller to generate training data
% Control gains for PD controller
Kp = 2.0;    % Proportional gain
Kd = 1.5;    % Derivative gain for reference trajectory
Kv = -0.5;   % Derivative gain for system velocity

% Generate training data using a PD controller
X = zeros(nx, length(t));
U = zeros(1, length(t));
X(:,1) = x0;
ref_velocity = zeros(1, length(t));

for k = 1:length(t)-1
    % PD control law components:
    error = ref_position(k) - X(1,k);    % Position error
    ref_velocity(k) = (ref_position(k+1)-ref_position(k))/dt;   % Reference velocity (numerical derivative)
    system_velocity = X(2,k);   % System velocity (state 2)

    % Full PD control law
    U(k) = Kp*error + Kd*ref_velocity(k) + Kv*system_velocity;

    % Optional: Add control input saturation
    U(k) = min(max(U(k), -2), 2);  % Limit control between -2 and 2

    % Simulate system one step forward
    X(:,k+1) = f(X(:,k), U(k));
end

% Visualization of PD control performance
figure('Name','Measurement Data')
subplot(3,1,1)
plot(t, ref_position, 'k--', t, X(1,:), 'b-')
legend('Reference', 'PD Response')
title('PD Controller Training Data Generation')
ylabel('Position')

subplot(3,1,2)
plot(t, X(2,:), 'b-')
title('System Velocity')
ylabel('Velocity')

subplot(3,1,3)
plot(t, U, 'r-')
title('Control Input')
ylabel('Control')
xlabel('Time')

%% Define specialized Koopman observables for control

% Observable functions library:
obsFcn = @(x,u,r) {
    x(1);                   % 1 State 1
    x(2);                   % 2 State 2
    r;                      % 3 Reference signal
    r - x(1);               % 4 Tracking error
    x(2)^2;                 % 5 Kinetic energy-like term
    u;                      % 6 Control input
    u^2;                    % 7 Control energy
    (r - x(1))^2;           % 8 Squared tracking error
    x(1)*x(2);              % 9 Mixed state term
    (r - x(1))*x(2);        % 10 Error-velocity product
    };

% Compute observables along trajectory
nobs = 10; %
Psi = zeros(nobs, length(t)); % initialize the obervable functions array
for k = 1:length(t)
    observed = obsFcn(X(:,k), U(k), ref_position(k));
    Psi(:,k) = cell2mat(observed);
end

%% Compute regularized Koopman operator

% Create data matrices with regularization
Psi_x = Psi(:,1:end-1);
Psi_y = Psi(:,2:end);
lambda = 1e-2; % regularization parameter
n = size(Psi_x,1);

% Compute regularized Koopman operator using ridge regression with Tikhonov
% regularization
% K = argmin ||Psi_y - K*Psi_x||^2 + lambda*||K||^2
K = Psi_y * Psi_x' / (Psi_x*Psi_x' + lambda*eye(n));

% Alternative formulation using SVD for better numerical stability
[svd_U,svd_S,svd_V] = svd(Psi_x, 'econ');
S_reg = diag(svd_S) ./ (diag(svd_S).^2 + lambda);
K_svd = Psi_y * (svd_V * diag(S_reg) * svd_U');

% Compute prediction errors for the two methods
error_ridge = norm(Psi_y - K*Psi_x, 'fro');
error_svd = norm(Psi_y - K_svd*Psi_x, 'fro');

% select best method
[~,idx] = min([error_ridge, error_svd]);
if idx == 2
    K = K_svd;
end

save("KoopmanMatrix.mat","K");

%% Compare Koopman predictions with measured data

% Preallocate predicted states
X_koopman = zeros(nx, length(t));  % nx = 2 if you consider just the actual system states
X_koopman(:,1) = X(:,1);           % Use the same initial condition as measured data

% For storing predicted observables
Psi_pred = zeros(nobs, length(t));
Psi_pred(:,1) = Psi(:,1);  % The first set of observables from training data

for k = 1:length(t)-1

    % Current predicted state (first 2 entries of Psi_pred)
    xk_pred = X_koopman(:,k);

    % Construct observables from predicted state, but use the same control and reference
    % as in training. (If you have a separate test set, you would replace U(k), ref_position(k), etc.)
    currentRef = ref_position(k);
    psi_k = obsFcn(xk_pred, U(k), currentRef);  % cell array from obsFcn
    psi_k = cell2mat(psi_k);

    % Predict next observables using Koopman operator
    psi_next = K * psi_k;

    % Store all predicted observables
    Psi_pred(:,k+1) = psi_next;

    % Extract the first 2 as the predicted states
    X_koopman(:,k+1) = psi_next(1:nx);
end

% Plot the comparison (e.g., compare x(1) and x(2))
figure('Name','Koopman Predictions vs Measured Data')
subplot(2,1,1)
plot(t, X(1,:), 'b-', t, X_koopman(1,:), 'r--', 'LineWidth',1.5)
legend('Measured x(1)','Koopman Pred x(1)','Location','Best')
title('Comparison of x(1) (Position)')

subplot(2,1,2)
plot(t, X(2,:), 'b-', t, X_koopman(2,:), 'r--', 'LineWidth',1.5)
legend('Measured x(2)','Koopman Pred x(2)','Location','Best')
title('Comparison of x(2) (Velocity)')
xlabel('Time (s)')

%% Design MPC using Koopman model with nlmpcMultistage

% Get the list of all installed toolboxes.
installedToolboxes = ver;

% Check if Model Predictive Control Toolbox is installed.
mpcInstalled = any(strcmpi('Model Predictive Control Toolbox', {installedToolboxes.Name}));

% Display results for Model Predictive Control Toolbox
if mpcInstalled
    disp('Model Predictive Control Toolbox is installed.');
    % Check if a license is available
    if license('test','MPC_Toolbox')
        disp('  A valid license for Model Predictive Control Toolbox is available.');
    else
        disp('  No valid license for Model Predictive Control Toolbox.');
    end
else
    disp('Model Predictive Control Toolbox is NOT installed.');
end

% Create nlmpcMultistage object
predHorizon = 10;  % Prediction horizon
nu = 1;  % Number of inputs

msobj = nlmpcMultistage(predHorizon, nx_extended, nu);

% Set sampling time
msobj.Ts = 0.01;

% Reference tracking target
ref = 1;  % reference position

% Define cost weights
Q = 10;   % Stage state cost weight
R = 0.1;  % Control cost weight
S = 20;   % Terminal state cost weight
Qi = 1;   % Integral state error weight

% Define state function parameters
pvstate = ref;

% Package parameters for cost function
pvcost = [ref; predHorizon; Q; R; S; Qi];

% Stage parameters (used by cost function)
pvstage = pvcost;

% Set control input constraints
msobj.ManipulatedVariables(1).Min = -5;
msobj.ManipulatedVariables(1).Max = 5;

% Define state and cost functions using MyNMPC class
clear MyNMPC
msobj.Model.StateFcn = "MyNMPC.StateFcn";

msobj.Model.ParameterLength = length(pvstate);  % pvstate length
msobj.Model.IsContinuousTime = false;

[msobj.Stages.CostFcn] = deal("MyNMPC.CostFcn");
[msobj.Stages.CostJacFcn] = deal("MyNMPC.CostJacFcn");
[msobj.Stages.ParameterLength] = deal(length(pvcost));

% Solver settings
msobj.Optimization.Solver = "fmincon";
msobj.Optimization.SolverOptions   = optimoptions('fmincon', 'Display', 'off');
msobj.Optimization.SolverOptions.Algorithm = "sqp";
msobj.Optimization.SolverOptions.MaxIterations = 20;
msobj.Optimization.SolverOptions.ConstraintTolerance = 1e-5;
msobj.Optimization.UseSuboptimalSolution = false;

% Check validity
simdata = getSimulationData(msobj);

simdata.StateFcnParameter = pvstate;
simdata.StageParameter = repmat(pvstage, predHorizon + 1, 1);
validateFcns(msobj,rand(nx_extended,1),rand(nu,1),simdata)

%% Simulate controlled system with Koopman MPC

% Define Reference Signal Parameters
ref_change_interval = 5;      % Interval in seconds to change reference
ref_min = 0.5;                 % Minimum reference value
ref_max = 1.5;                 % Maximum reference value

% Simulation Parameters
Tf = 20;                      % Total simulation time in seconds
dt = msobj.Ts;                % Sampling time, e.g., 0.01
t = 0:dt:Tf - dt;             % Time vector

numSteps = length(t)-1;       % Number of simulation steps

% Generate Step-Like Reference Signal
rng("default") % for reproducible results
num_changes = ceil(Tf / ref_change_interval);
ref_values = ref_min + (ref_max - ref_min) * rand(num_changes, 1);  % Random references

% Initialize reference vector
ref_signal = zeros(1, length(t));

for i = 1:num_changes
    start_time = (i-1)*ref_change_interval;
    end_time = min(i*ref_change_interval, Tf);

    % Find indices corresponding to the current interval
    idx = find(t >= start_time & t < end_time);

    % Assign the reference value for this interval
    ref_signal(idx) = ref_values(i);
end

% Handle any remaining time if Tf is not a multiple of ref_change_interval
if i*ref_change_interval < Tf
    idx = find(t >= i*ref_change_interval);
    ref_signal(idx) = ref_min + (ref_max - ref_min) * rand(1);
end

%% MATLAB simulation loop

% Initial conditions
x0 = [0; 0; 0]; % [x1; x2; xi] Original states and integral error

% Initialize storage
X_mpc = zeros(nx_extended, length(t));
U_mpc = zeros(1, length(t));
X_mpc(:,1) = x0;
iter = zeros(1, length(t));

% Initialize MPC controller simdata struct
simdata = getSimulationData(msobj);

% First step control, u(0)
uk = 0;

disp("Simulate with MATLAB...")

for k = 1:numSteps
    % Current state vector
    xk = X_mpc(:,k);

    % Current reference
    current_ref = ref_signal(k);

    % Update simdata
    simdata.StateFcnParameter = current_ref;
    pvcost = [current_ref; predHorizon; Q; R; S; Qi];
    simdata.StageParameter = repmat(pvcost, predHorizon + 1, 1);

    % Compute optimal control move
    [uk,simdata,info] = nlmpcmove(msobj, xk, uk, simdata);

    % Solver iterations
    iter(k) = info.Iterations;

    % Apply control and simulate system
    U_mpc(k) = uk;

    % Simulate system using the "true" model for [x1; x2]
    x_next_physical = f(X_mpc(1:2,k), U_mpc(k));   % [x1_next; x2_next]

    % Compute tracking error
    e_k = current_ref - X_mpc(1,k);  % e(k) = ref - x1(k)

    % Update integral of the error
    xI_next = X_mpc(3,k) + e_k * dt;  % xi(k+1) = xi(k) + e(k)*dt

    % Combine physical states and integral state
    X_mpc(:,k+1) = [x_next_physical; xI_next];

    % Display progress every 100 steps
    if mod(k,100) == 0
        disp("Step " + k + "/" + numSteps + ", Cost " + info.Cost + ", Position " + X_mpc(1,k+1))
    end
end

%% Plot results
figure('Name','Control performance')
subplot(5,1,1)
hold on
plot(t, ref_signal, 'k--', 'LineWidth', 1.5)           % Reference signal
plot(t, X_mpc(1,:), 'b-', 'LineWidth',1.5)            % Controlled position
legend('Reference', 'Controlled Output')
title('Tracking Performance')
ylabel('Position')
ylim([min(ref_min, min(X_mpc(1,:))) - 0.1, max(ref_max, max(X_mpc(1,:))) + 0.1])

subplot(5,1,2)
plot(t, X_mpc(2,:), 'b-', 'LineWidth',1.5)
title('Velocity')
ylabel('Velocity')

subplot(5,1,3)
plot(t, X_mpc(3,:), 'g-', 'LineWidth',1.5)
title('Integral of Tracking Error')
ylabel('Integral Error')

subplot(5,1,4)
plot(t, U_mpc, 'r-', 'LineWidth',1.5)
title('Control Input')
ylabel('Control')

subplot(5,1,5)
plot(t, iter, 'r-', 'LineWidth',1.5)
title('Solver iterations')
xlabel('Time (s)')
ylabel('Iterations')

%% Using Simulink

% To speed up computation time during simulation, use getCodeGenerationData
% and buildMEX to generate code for the Multistage MPC object with the
% specified properties. (MATLAB Coder is needed) 

% Check if MATLAB Coder is installed.
matlabCoderInstalled = any(strcmpi('MATLAB Coder', {installedToolboxes.Name}));

% Display results for MATLAB Coder
if matlabCoderInstalled
    disp('MATLAB Coder is installed.');
    % Check if a license is available
    if license('test','MATLAB_Coder')
        disp('  A valid license for MATLAB Coder is available.');
    else
        disp('  No valid license for MATLAB Coder.');
    end
else
    disp('MATLAB Coder is NOT installed.');
end


% Initial conditions
x0 = [0; 0; 0]; % [x1; x2; xi] Original states and integral error
u0 = zeros(nu,1);

% Generate data structures using getCodeGenerationData.
[coredata,onlinedata] = getCodeGenerationData(msobj,x0,u0);

% Generate a MEX function called nlmpcControllerMEX.
MEXname = 'nlmpcControllerMEX_withIntegralControl';
mexFcn = buildMEX(msobj, MEXname, coredata,onlinedata);
 
% open controller model
load_system("multistageNLMPC_controlBlock_withIntegralControl_mdlRef");

% set MEX file
set_param("multistageNLMPC_controlBlock_withIntegralControl_mdlRef/Multistage Nonlinear MPC",mexname = MEXname);

% open simulation model
open_system("test_nlmpcMultistage_withIntegralControl");

% open scope
open_system("test_nlmpcMultistage_withIntegralControl/Scope")

disp("Simulate with Simulink...")

% simulate
out = sim("test_nlmpcMultistage_withIntegralControl");

% check simulation speed using the ratio: simulation time/elapsed wall time
simSpeed = (out.SimulationMetadata.ModelInfo.StopTime - out.SimulationMetadata.ModelInfo.StartTime)/...
out.SimulationMetadata.TimingInfo.ExecutionElapsedWallTime ;

disp("Simulation speed using the ratio: simulation time/elapsed wall time: "+simSpeed)