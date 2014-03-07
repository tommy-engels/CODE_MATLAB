function dt = dt_EPS( eps )
%% dt = dt_EPS( eps ) returns the time step allowed by explicit penalization
% requires params.CFL, params.dx
    global params
    dt = 0.90*eps;
end