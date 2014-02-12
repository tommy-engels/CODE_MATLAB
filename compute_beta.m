function beta = compute_beta (u)
    global params
    ux_x = params.D1x*u(:,:,1) ;
    ux_y = u(:,:,1)*params.D1y' ;
    uy_x = params.D1x*u(:,:,2) ;
    uy_y = u(:,:,2)*params.D1y' ;
    
    % actual beta field (normal derivatives, component-wise)
    beta(:,:,1) = (params.n_x.*ux_x + params.n_y.*ux_y);
    beta(:,:,2) = (params.n_x.*uy_x + params.n_y.*uy_y);
end


%     ux_x = cofitxy(1i*params.Kx.*fft2(u(:,:,1)));
%     ux_y = cofitxy(1i*params.Ky.*fft2(u(:,:,1)));
%     uy_x = cofitxy(1i*params.Kx.*fft2(u(:,:,2)));
%     uy_y = cofitxy(1i*params.Ky.*fft2(u(:,:,2)));
    
%     [ux_x,ux_y,uy_x,uy_y] = upwind_differences(u);