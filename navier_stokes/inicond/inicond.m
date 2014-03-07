function [u, uk] = inicond()
%% [u0, uk] = inicond() initial conditions, according to the flag params.inicond
global params

switch params.inicond
    case 'quiescent'
        u = zeros(params.nx,params.ny,2);
        uk = coftxy_2d(u);
end

fprintf('initial field divergence: %12.5e\n',max(max(abs(cofitxy(divergence_2d(uk))))));
end
