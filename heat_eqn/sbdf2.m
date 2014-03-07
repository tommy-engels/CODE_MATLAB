function [u_new, uk_new] = sbdf2(time,it,u_old,uk_old,u_oldold,uk_oldold)
    %% SBDF2 scheme
    % note here one term is explicit (the diffusion) and one implicit (the
    % penalization). background is that you hoped to have a smaller
    % splitting error using SBDF2, because you'd treat the nonlinear term
    % explicitly anyways. Didn't work; shows the same behavior
    global params
    dt = params.dt;
    nu = params.nu;
    
    if (it==1)
        [u_new, uk_new] = strang_spectral(time,it,u_old,uk_old,u_oldold,uk_oldold);
    else    
    % laplace operator
    K2 = -params.Kx.^2 - params.Ky.^2;
    
    uk_old = dealias(uk_old);
    uk_oldold = dealias(uk_oldold);
    
    % right hand side    
    B = cofitxy(2*dt*nu*K2.*uk_old - dt*nu*K2.*uk_oldold) + 2*u_old - 0.5*u_oldold;
    
    % "solve" the diagonal system (that's really easy!)
    u_new = B ./ ( 3/2 + params.dt*params.mask/params.eta);
    
    % consistent output
    uk_new = dealias(fft2( u_new ));
    u_new = cofitxy( uk_new );
    end
end