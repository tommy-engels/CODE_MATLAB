function [u0, uk] = inicond()
%% initial condition for lamballais case
    global params
     
          
    u0 = params.u_ex*0;
    uk = coftxy_2d(u0);
    uk = make_incompressible_2d (uk);
    u0 = cofitxy_2d(uk);

    
    % go to Fourier space
    uk = coftxy_2d(u0);
    fprintf('initial field divergence: %12.5e\n',max(max(abs(cofitxy(divergence_2d(uk))))));
end