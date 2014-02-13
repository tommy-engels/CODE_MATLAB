function [n_x, n_y] = compute_normals( phi )
    %% computes the normal vector field from the signed distance funciton
    global params
    % differentiate phi field numerically    
    phi_dx = cofdx_fd ( phi, D1(params.nx,params.dx) );
    phi_dy = cofdy_fd ( phi, D1(params.ny,params.dy) );
    
    % normalization
    normgrad = sqrt( phi_dx(:,:).^2 + phi_dy(:,:).^2  );
    normgrad ( abs(normgrad) < 1e-10 ) = 1.0;
    
    % these are the normal vectors:
    n_x = phi_dx(:,:) ./ normgrad;
    n_y = phi_dy(:,:) ./ normgrad;        
    
    
    %% when using chantalat: reduce normals to a boundary layer
    % compute boundary-layer-blending function (only when doing
    % chantalat-type prolongation)
    if strcmp(params.CASE,'chantalat')
        params.blend = zeros(params.nx,params.ny);
        for ix=1:params.nx
            for iy=1:params.ny
                params.blend(ix,iy) = smoothstep( abs(phi(ix,iy)), 4*params.dx, 3*params.dx);
            end
        end
        
        % blend normals to a boundary layer in the vicinity of the
        % fluid/solid interface
        n_x = n_x.*params.blend;
        n_y = n_y.*params.blend;
    end    
end