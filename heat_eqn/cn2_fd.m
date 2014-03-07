function [u_new, uk_new] = cn2_fd(time,it,u,uk,u_oldold,uk_oldold)
    %% ordinary CN2, treating all terms implicitly, no splitting, FD only
    % all terms are treated with the same implicit integrator. This scheme
    % allows eta<<dt and indeed does not show saturation: the color plot is
    % almost vertical lines. this means for a given dt, you CAN decrease
    % eta and really compute the error accordingly. no saturation.
    global params
    chi = params.mask/params.eta;
        
    D = sparse(D24p(params.nx,params.dx));
    if ( it == 1 ) || ( params.T_end-time < params.dt ) % do this for the last step as well
        I = speye(params.nx);
        A = kron(D,I) + kron(I,D);
        I = speye( (params.nx)*(params.ny) );
        M = I - 0.5*params.dt*(params.nu*A - sparse(diag(reshape(chi,[],1))) );
        [params.L, params.U, params.P, params.Q, params.R] = lu( M );        
    end
    
    b = u + 0.5*params.dt*( params.nu*(cofdx_fd(u,D)+cofdy_fd(u,D)) - chi.*u);
    
    u_new = params.Q * (params.U \ (params.L \ (params.P * (params.R \ reshape( b ,[],1))))) ;
    u_new = reshape(u_new,params.nx,params.ny);
    
    uk_new = u_new; % dummy argument; discarded
end