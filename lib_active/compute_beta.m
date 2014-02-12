function beta = compute_beta (u)
    %% compute BETA field
    % the beta field is the "normal" derivative of the velocity field. so its
    % actually d/dn of each component of u. The n-direction is indicated by the
    % normal vectors, defined as the gradient of phi.  
    global params
    
    if strcmp(params.active_beta,'central')
        ux_x = params.D1x*u(:,:,1) ;
        ux_y = u(:,:,1)*params.D1y' ;
        uy_x = params.D1x*u(:,:,2) ;
        uy_y = u(:,:,2)*params.D1y' ;
        
    elseif strcmp(params.active_beta,'upwind')
        [ux_x,ux_y,uy_x,uy_y] = upwind_differences(u);
        
    elseif strcmp(params.active_beta,'spectral')        
        ux_x = cofitxy(1i*params.Kx.*fft2(u(:,:,1)));
        ux_y = cofitxy(1i*params.Ky.*fft2(u(:,:,1)));
        uy_x = cofitxy(1i*params.Kx.*fft2(u(:,:,2)));
        uy_y = cofitxy(1i*params.Ky.*fft2(u(:,:,2)));
    end
    
    % actual beta field (normal derivatives, component-wise)
    beta(:,:,1) = (params.n_x.*ux_x + params.n_y.*ux_y);
    beta(:,:,2) = (params.n_x.*uy_x + params.n_y.*uy_y);
end



function  [ux_x,ux_y,uy_x,uy_y] = upwind_differences(u)
    % chantalat uses upwind differences to compute the normal derivative field,
    % we shall try to do so as well.
    global params
    dx = params.dx;
    dy = params.dy;
    
    ux_x = zeros(params.nx,params.ny);
    ux_y = zeros(params.nx,params.ny);
    uy_x = zeros(params.nx,params.ny);
    uy_y = zeros(params.nx,params.ny);
    
    for ix=2:params.nx-1
        for iy=2:params.ny-1
            % derivative in x-direction
            if (params.phi(ix+1,iy)<params.phi(ix-1,iy))
                ux_x(ix,iy) = (u(ix+1,iy,1)-u(ix,iy,1)) / dx;
                uy_x(ix,iy) = (u(ix+1,iy,2)-u(ix,iy,2)) / dx;
            else
                ux_x(ix,iy) = (u(ix,iy,1)-u(ix-1,iy,1)) / dx;
                uy_x(ix,iy) = (u(ix,iy,2)-u(ix-1,iy,2)) / dx;
            end
            % derivative in y-direction
            if (params.phi(ix,iy+1)<params.phi(ix,iy-1))
                ux_y(ix,iy) = (u(ix,iy+1,1)-u(ix,iy,1)) / dy;
                uy_y(ix,iy) = (u(ix,iy+1,2)-u(ix,iy,2)) / dy;
            else
                ux_y(ix,iy) = (u(ix,iy,1)-u(ix,iy-1,1)) / dy;
                uy_y(ix,iy) = (u(ix,iy,2)-u(ix,iy-1,2)) / dy;
            end
        end
    end
end