%% irf_func - calculates IRF
%
% This function uses the "Dynare" steady state (i.e. variables are in
% deviations from non-linear steady state).
% =========================================================================
function out = irf_func(oo_, mats, e1, horizon, bypass_Callum) 

Q = mats.Q ; 
G = mats.G ; 

if bypass_Callum
   x = [oo_.steady_state] ;
else 
    x = [oo_.steady_state;0] ;
end

y(:,1) = Q * x + G * e1 ;

for t_ = 2 : horizon 
    x       = y(:,t_-1) ; 
    y(:,t_) = Q * x ;
end

out = y ;
