%% tv_mats_func.m
%
% This function calculates time varying reduced form solution matrices
% using the approach described in:
% Kulish and Pagan (2017), "Estimation and Solution of Models with
% Expectations and Structural Changes", Journal of Applied Econometrics
% 32(2), pp. 255-274.
%
% For notation purposes, the linearized model is:
% A x(t) = B x(t-1) + D E_t x(t+1) + E e(t) 
% 
% =========================================================================
function out = tv_mats_func(input, mats_init, mats_final, M_) ;

% Inputs 
horizon     = input.horizon ;   % Length of the simulation
t_a         = input.t_a ;       % Date of structural change
t_b         = input.t_b ;       % Date of anticipation / realisation of structural change

n_bp        = M_.endo_nbr + 1 ;     % Number of endogenous variables (+ 1 because of constant)
n_shock     = M_.exo_nbr ;          % Number of exogenous variables

% Initialise matrices
Q_tv = zeros(n_bp, n_bp, horizon) ;
G_tv = zeros(n_bp, n_shock, horizon) ;
Q_init = mats_init.Q ;

% =========================================================================
% Three cases
% =========================================================================

% ** Case 1 = date of structural change = date of realisation **
if t_a == t_b

Q_tv(:,:,1 : t_a-1) = repmat(mats_init.Q,1,1,t_a-1) ;
G_tv(:,:,1 : t_a-1) = repmat(mats_init.G,1,1,t_a-1) ;

Q_tv(:,:, t_a : horizon) = repmat(mats_final.Q,1,1,horizon - t_a + 1) ;
G_tv(:,:, t_a : horizon) = repmat(mats_final.G,1,1,horizon - t_a + 1) ;

% ** Case 2 => Anticipated structural change **
elseif t_b < t_a

% Before announcement of structural change    
Q_tv(:,:,1 : t_b-1) = repmat(mats_init.Q,1,1,t_b-1) ;
G_tv(:,:,1 : t_b-1) = repmat(mats_init.G,1,1,t_b-1) ;

% After implementation of structural change
Q_tv(:,:, t_a : horizon) = repmat(mats_final.Q,1,1,horizon - t_a + 1) ;
G_tv(:,:, t_a : horizon) = repmat(mats_final.G,1,1,horizon - t_a + 1) ;

% Intermediate periods => Initial structural matrices but anticipation of
% future changes
Bt = mats_init.A \ mats_init.B ;
Dt = mats_init.A \ mats_init.D ;
Et = mats_init.A \ mats_init.E ;

for t_ = t_a-1 : -1 : t_b
Q_tv(:,:,t_) = (eye(n_bp) - Dt * Q_tv(:,:,t_+1)) \ Bt ;
G_tv(:,:,t_) = (eye(n_bp) - Dt * Q_tv(:,:,t_+1)) \ Et ;
end

% ** Case 3 =>  Date of structural change before realisation **
elseif t_a < t_b

% Before structural change
Q_tv(:,:,1 : t_a-1) = repmat(mats_init.Q,1,1,t_a-1) ;
G_tv(:,:,1 : t_a-1) = repmat(mats_init.G,1,1,t_a-1) ;

% After realisation of structural change
Q_tv(:,:, t_b : horizon) = repmat(mats_final.Q,1,1,horizon - t_b + 1) ;
G_tv(:,:, t_b : horizon) = repmat(mats_final.G,1,1,horizon - t_b + 1) ;

% Intermediate periods => Final structural matrices, expectations based on
% initial matrices
Bt = mats_final.A \ mats_final.B ;
Dt = mats_final.A \ mats_final.D ;
Et = mats_final.A \ mats_final.E ;

for t_ = t_b-1 : -1 : t_a
Q_tv(:,:,t_) = (eye(n_bp) - Dt * Q_init) \ Bt ;
G_tv(:,:,t_) = (eye(n_bp) - Dt * Q_init) \ Et ;
end

end

% Outputs
out.Q_tv = Q_tv;
out.G_tv = G_tv ;

