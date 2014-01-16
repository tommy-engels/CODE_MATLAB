function ExactSolution_Couette
    global params    
    x0 = 0.5*params.Lx;
    y0 = 0.5*params.Ly;
    
    create_mask();
        
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % EXACT SOLUTION
    %%%%%%%%%%%%%%%%%%%%%%%%%
    R = sqrt( (params.X-x0).^2 + (params.Y-y0).^2);
    
    % abs velocity exact
    A = -params.omega*params.R1^2 / (params.R2^2 - params.R1^2);
    B = +params.omega*(params.R1^2)*(params.R2^2) / (params.R2^2 - params.R1^2);    
    u_ex(:,:) = (A*R + B./R) .* (1-params.mask);
    
    % radial unity vector
    er(:,:,2) = +(params.X-x0)./R;
    er(:,:,1) = -(params.Y-y0)./R;
    
    % exact solution
    params.u_ex = zeros(params.nx,params.ny,2);
    params.u_ex(:,:,1) = er(:,:,1).*u_ex ;
    params.u_ex(:,:,2) = er(:,:,2).*u_ex ;
    
    % delete velocity inside body
    for ix=1:params.nx
        for iy=1:params.ny
        if (params.mask(ix,iy)>0)
            params.u_ex(ix,iy,1:2) = 0;
        end
        end
    end
    
    % set solid body rotation in the smaller cylinder (used in classic
    % penalization)
    for ix=1:params.nx
        for iy=1:params.ny
        if (R(ix,iy)<=params.R1)
            params.u_ex(ix,iy,1) =-params.omega*(params.Y(ix,iy)-y0);
            params.u_ex(ix,iy,2) = params.omega*(params.X(ix,iy)-x0);
        end
        end
    end
 end