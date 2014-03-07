function [u_new, uk_new] = strang_fd(time,it,u,uk,u_oldold,uk_oldold)
    %% strang splitting with finite differences for the laplacian + CN2
    global params
    chi = params.mask; us=params.us;
    %----------------------------------------------------------------------
    % 1st strang step, half a time step for penalization term
    %----------------------------------------------------------------------
    dt1 = 0.5*params.dt; 
    a = dt1/params.eta;    
    u_new  = (u-chi.*us).*exp(-chi*a) + chi.*us;    
    %----------------------------------------------------------------------
    % 2nd strang step, RK2 time step for NS
    %----------------------------------------------------------------------
    D = sparse(D24p(params.nx,params.dx));
    if ( it == 1 ) || ( params.T_end-time < params.dt ) % do this for the last step as well
        I = speye(params.nx);
        A = kron(D,I) + kron(I,D);
        I = speye( (params.nx)*(params.ny) );
        M = I - params.nu*0.5*params.dt*A;
        [params.L, params.U, params.P, params.Q, params.R] = lu( M );        
    end
    b = u_new + 0.5*params.nu*params.dt*( cofdx_fd(u_new,D)+cofdy_fd(u_new,D) );
    
    u_new = params.Q * (params.U \ (params.L \ (params.P * (params.R \ reshape( b ,[],1))))) ;
    u_new = reshape(u_new,params.nx,params.ny);    
    %----------------------------------------------------------------------
    % 3rd strang step, half a time step for penalization term
    %----------------------------------------------------------------------
    dt1 = 0.5*params.dt; 
    a = dt1/params.eta;  
    u_new  = (u_new-chi.*us).*exp(-chi*a) + chi.*us;
    uk_new = u_new; % dummy argument; discarded
end