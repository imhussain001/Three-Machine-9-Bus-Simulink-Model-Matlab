%WSCC9_LFC_INIT  Workspace for wscc9_lfc_governor.slx (load-frequency study).
%   WSCC 3-machine 9-bus system, classical SGs (model 0.0) as in
%   wscc9_transient_stability.slx, each with a simple non-reheat governor:
%   dPm = -(1/R) dw -> 1/(1+s Tg) -> 1/(1+s Tt). The load step is a
%   constant-impedance conductance added at one load bus, i.e. the
%   Kron-reduced network switches from Ypre to Ypost at t_step.
%   Runs automatically as the model's InitFcn. Powers on the 100 MVA base.

% Network, internal voltages E1..E3, initial angles D10..D30, H, f
evalc('wscc9_transient_init');

%% Load step
dP_L_MW  = 20;     % [MW] load change, + increase / - rejection
step_bus = 5;      % load bus: 5, 6 or 8
t_step   = 1;      % [s]

load_buses = [5 6 8];
V_load     = [0.996 1.013 1.016];   % load-flow voltages used for the loads in wscc9_transient_init
assert(ismember(step_bus, load_buses), 'step_bus must be one of the load buses 5, 6, 8');
S_base = 100;      % [MVA]
dG_L   = (dP_L_MW/S_base) / V_load(load_buses == step_bus)^2;   % conductance at load-flow voltage

Ystep = Ybusext;
Ystep(step_bus, step_bus) = Ystep(step_bus, step_bus) + dG_L;
Ypre  = Y1;                              % intact network, nominal load
Ypost = kron_reduce(Ystep, 3);           % generator internal nodes 1..3 kept

%% Operating point of the model
Eph = [E1; E2; E3] .* exp(1j*[D10; D20; D30]);
Pe0 = real(Eph .* conj(Ypre*Eph));       % [pu, 100 MVA] electrical power at t = 0
Pm0 = Pe0;                               % equilibrium: Pm = Pe (absorbs load-flow rounding)
fprintf('Pm0 = [%s] MW (load-flow Pg = [71.6 163 85] MW)\n', num2str(100*Pm0', '%8.2f'));

%% Simple governor (non-reheat), droop on own machine base
S_MVA = [247.5; 192; 128];               % machine ratings (H above is already on 100 MVA)
kS    = S_MVA / S_base;
R     = 0.05;                            % droop [pu own]
Tg    = 0.2;                             % governor time constant [s]
Tt    = 0.5;                             % turbine time constant [s]
ws    = 2*pi*f;                          % [rad/s]

% Damping in the swing equation, D*dw on own base. The classical SG has no
% damping, so its 1.40 Hz and 2.14 Hz inter-machine modes are undamped; with
% the governors they become unstable (zeta = -1.2 %, -0.7 % at D = 0).
% D = 4 gives zeta = +1.7 %, +2.3 %. D also adds to the frequency response.
D     = 4;

%% Expected steady state (primary control only, frequency-independent loads)
beta    = sum(kS) * (1/R + D);           % [pu_100 / pu_f], droop + damping
df_pred = -(dP_L_MW/S_base) / beta;      % [pu]
% The constant-impedance load draws less than dP_L_MW as its voltage drops,
% so the simulated deviation is somewhat smaller than this nominal estimate.
fprintf('Nominal steady-state frequency for %+g MW: %.4f Hz (df = %.4f Hz)\n', dP_L_MW, f*(1 + df_pred), f*df_pred);

function Yr = kron_reduce(Y, n)
%KRON_REDUCE  Eliminate all nodes after the first n.
K = Y(1:n, 1:n);  L = Y(1:n, n+1:end);  M = Y(n+1:end, n+1:end);
Yr = K - L*(M\L.');
end
