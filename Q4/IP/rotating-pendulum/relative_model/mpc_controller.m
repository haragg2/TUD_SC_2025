close all
addpath('functions\');
addpath('linear_system_matrices\');

%% Linearize the system

sys_mat = load("LinearSys_pi_0.mat");
sys_ct = ss(sys_mat.sys_A, sys_mat.sys_B, [1, 0, 0, 0; 0, 0, 1, 0], []);

%% Discretize the system

sys_dis = c2d(sys_ct, h);

% System matrices
A = sys_dis.A;
B = sys_dis.B;
C = sys_dis.C;
D = sys_dis.D;

%% Bounds

xlb = [-1.5; -15; -pi; -15];
xub = -xlb;
ulb = -1;
uub = 1;

N = 12; % Horizon

% Defines the dimensions
dim.nx = size(sys_dis.A, 1);        % Number of states
dim.nu = size(sys_dis.B, 2);        % Number of inputs
dim.ny = size(sys_dis.C, 1);        % Number of outputs

%% LQR

% Tuning weights

R = 0.1;

Q = [40 0 0 0;
     0 10 0 0;
     0 0 70 0;
     0 0 0 10];

[K_lqr, P, ~] = dlqr(A, B, Q, R);
K_lqr = [K_lqr, 0.02, -0.02];

K = -K_lqr(1:4);

% Controllability check
ctrb_sys = ctrb(A, B);
unco = length(A) - rank(ctrb_sys);
if unco > 0
    warning('Discretized linear model is uncontrollable');
end

% Observability check
obs_sys = obsv(A, C);
unob = length(A) - rank(obs_sys);
if unob > 0
    warning('Discretized linear model is unobservable');
end

% Check if (A, Q') has no unobservable modes on the unit circle
isNoUnobservableModes = rank(ctrb(A, transpose(Q))) == size(A,1);
if isNoUnobservableModes == 0
    warning("Pair (A, Q') has some unobservable modes");
end

%% Compute X_f

Xn = struct();
V = struct();
Z = struct();

[Xn.('lqr'), V.('lqr'), Z.('lqr')] = findXn(A, B, K, N, xlb, xub, ulb, uub, 'lqr');
model_mpc = struct('A', A, 'B', B, 'C', C, 'Bd', zeros(size(B)), 'Cd', zeros(size(C, 1), 1), 'N', N);
constraint = Z.lqr;
penalty = struct('Q', Q, 'R', R, 'P', P);
terminal = Xn.lqr{1}; % LQR terminal set

%% Regulation MPC

xr = [0; 0; 0; 0];
model_matrices = buildmatrices(xr, model_mpc, constraint, penalty, terminal);

%% Reference Tracking - Comment the section out if running Regulation MPC

x_ref = pi/6;
yref = [-x_ref; x_ref];
d_hat = 0;

xr = targetSelector(model_mpc, Z.lqr, dim, d_hat, yref);
model_matrices = buildmatrices(xr, model_mpc, constraint, penalty, terminal);


%% Functions

function xr = targetSelector(LTI, Z, dim, d_hat, yref)

    eqconstraints.A = [eye(dim.nx) - LTI.A, -LTI.B; LTI.C, zeros(size(LTI.C, 1), dim.nu)];
    eqconstraints.b = [LTI.Bd * d_hat; yref - (LTI.Cd*d_hat)];

    ineqconstraints.A = [Z.('G'), Z.('H')];
    ineqconstraints.b = Z.('psi');

    H = blkdiag(zeros(dim.nx), eye(dim.nu));
    h = zeros(dim.nx+dim.nu, 1);

    options1 = optimoptions(@quadprog);
    options1.OptimalityTolerance=1e-20;
    options1.ConstraintTolerance=1.0000e-15;
    options1.Display='off';
    xur=quadprog(H,h,ineqconstraints.A,ineqconstraints.b,eqconstraints.A,eqconstraints.b,[],[],[],options1);
    xr = xur(1:dim.nx);
%     ur = xur(dim.nx+1:end);
end

function matrices = buildmatrices(xr,model, constraint, penalty, terminal)

% Returns the matrices struct used by the main function

N = model.N;
A = model.A;

if diff(size(A)) ~= 0
    error('model.Ak must be a square matrix.');
end
numx = size(A,1); % Number of states

B = model.B;
numu = size(B,2); % Number of inputs

if isfield(constraint,'G')
    G = constraint.G;
    H = constraint.H;
    psi = constraint.psi;
else
    G = zeros(0, numx);

    H = zeros(0, numu);
    psi = zeros(0, 1);
end

% Penalty matrices
Q = penalty.Q;
R = penalty.R;
P = penalty.P;
M = zeros(numx, numu);

littleH = [Q, M; M', R];
bigH = kron(eye(N),littleH);
bigH = blkdiag(bigH, P); % Add final penalty matrix

bigf = zeros(size(bigH, 1), 1);

% Structure of big A (both Aeq and Alt) is
% +-                  -+
% |  A1 A2  0  0  ...  |
% |   0 A1 A2  0  ...  |
% |   0  0 A1 A2  ...  |
% |  ...         A2  0 |
% |  ...         A1 A2 |
% +-                  -+
%
% We construct it first as
%
% +-             -+
% |  ... A1 A2  0 |
% |  ...  0 A1 A2 |
% |  ...  0  0 A1 | <= Note exta A1 that has to be removed.
% +-             -+
%
% and then get rid of extra rows (for last A1) and columns (for
% nonexistant variable u_N)

% For Equalities:
% A1 = [A, B], A2 = [-I, 0]
littleAeq1 = [A, B];
littleAeq2 = [-eye(size(A)), zeros(size(B))];

bigAeq = kron(eye(N+1),littleAeq1) + kron(diag(ones(1,N),1),littleAeq2);

bigAeq = bigAeq(1:end-numx,1:end-numu);
    % Remove columns for u_N and for rows of x_N+1 = A x_N + B u_N

bigbeq = zeros(numx*N,1);

% For Inequalities:
% A1 = [G, H], A2 = [0, 0]
littleAlt1 = [G, H];
littleAlt2 = [zeros(size(G)), zeros(size(H))];

bigAlt = kron(eye(N+1),littleAlt1) + kron(diag(ones(1,N),1), ...
    littleAlt2);
bigAlt = bigAlt(1:end-size(littleAlt1,1),1:end-numu);
    % Remove columns for u_N and for rows of G x_N + D u_N <= d
bigblt = repmat(psi, N, 1);

% Variable bounds
numVar = length(bigf);
LB = -inf*ones(numVar,1);
UB = inf*ones(numVar,1);

% Decide terminal constraint
if isstruct(terminal) && all(isfield(terminal, {'A', 'b'}))
    Af = terminal.A;
    bf = terminal.b + Af*xr;
    bigAlt = [bigAlt; [zeros(size(Af, 1), numVar - numx), Af]];
    bigblt = [bigblt; bf];
elseif isvector(terminal) && length(terminal) == numx
    LB(end-numx+1:end) = terminal;
    UB(end-numx+1:end) = terminal;
elseif isempty(terminal)
    % No terminal constraint. Pass
else
    error('Unknown input for terminal!');
end

% Build struct with appropriate names
matrices = struct('H', bigH, 'f', bigf, 'Aeq', bigAeq, 'beq', bigbeq, ...
                  'Alt', bigAlt, 'blt', bigblt, 'lb', LB, 'ub', UB);

end

