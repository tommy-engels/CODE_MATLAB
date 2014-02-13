function create_mask
    %% create mask (lamballais case)
    % ---------------------------------------------------------------------
    % this version is based on level set: it creates not only the mask, but
    % also the signed distance function phi and the normal vectors.
    % TASKS
    % * create the signed distance function phi
    % * create the mask function chi
    % * compute the normal vectors (possibly blend them to a small boundary
    %   layer in the vicinity of the interface)
    % * create the sponge mask (if that is used)
    % ---------------------------------------------------------------------
    global params
    % allocate memory
    params.mask = zeros(params.nx, params.ny);
    params.us   = zeros(params.nx,params.ny,2);
    params.phi = zeros(params.nx,params.ny);    
    params.masksponge = zeros(params.nx,params.ny);
    
    x0 = 0.5*params.Lx;
    y0 = 0.5*params.Ly;
    
    R = sqrt( (params.X-x0).^2 + (params.Y-y0).^2);
    params.R = R;
    
    load_lamballais_us()
    fprintf('grid points on the cylinder: %i\n',floor(params.R1*2/params.dx))
    
    %% compute signed distance function
    phi1 =  R - params.R1;
    phi2 =-(R - params.R2);
    phi3 = (R - params.R3);    
    
    for ix=1:params.nx
        for iy=1:params.ny
            v=[phi1(ix,iy) phi2(ix,iy) phi3(ix,iy)];
            [C,i]=min(abs(v));
            params.phi(ix,iy)=v(i);
        end
    end
    
    %% heaviside function (chi function)
    params.mask (params.phi <= 0.0) = 1.0;    

    %% sponge mask    
    if (strcmp(params.sponge,'yes'))
        params.masksponge ( (params.R>=params.R3+params.dx*5)&(params.X<=params.Lx/2) ) = 10.0;
    end

    %% compute normals
    [params.nx, params.ny] = compute_normals( params.phi );

   
    %% construct u_BC used to enforce non-homogeneous boundary conditions
    if strcmp(params.CASE,'chantalat')
        params.u_BC = params.us;
        kill = ones(size(params.mask));
        kill(R<=1.25*params.R1) = 0;
        params.u_BC(:,:,1) = params.u_BC(:,:,1) .* kill;
        params.u_BC(:,:,2) = params.u_BC(:,:,2) .* kill;
    else
        kill = ones(size(params.mask));
        kill(R<=1.25*params.R1) = 0;
        params.u_BC(:,:,1) = params.u_ex(:,:,1) .* kill;
        params.u_BC(:,:,2) = params.u_ex(:,:,2) .* kill;
    end
end