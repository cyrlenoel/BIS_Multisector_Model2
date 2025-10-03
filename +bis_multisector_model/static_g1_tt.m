function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 352);

T = bis_multisector_model.static_resid_tt(T, y, x, params);

T(345) = (-params(1297));
T(346) = (-(1/params(1340)*params(1345)*params(1335)));
T(347) = (-(T(2)*(params(1307)-T(123))));
T(348) = (-(1/params(1351)));
T(349) = 1-(1-params(1308));
T(350) = (-params(1299));
T(351) = (-params(1301));
T(352) = (-(1/params(1309)));

end
