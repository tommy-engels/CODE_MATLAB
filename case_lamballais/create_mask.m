function create_mask
    global params

    % allocate
    params.mask = zeros(params.nx, params.ny);
    params.us   = zeros(params.nx,params.ny,2);
    
    x0 = 0.5*params.Lx;
    y0 = 0.5*params.Ly;
    X = params.X;
    Y = params.Y;  
    R = sqrt( (params.X-x0).^2 + (params.Y-y0).^2);       

    %-----------------------------------
    % LAMBAILLAIS flow past cylinder
    %-----------------------------------
    load_lamballais_us()
    fprintf('grid points on the cylinder: %i\n',floor(params.R1*2/params.dx))
    
    % set the actual cylinder
    params.mask ( R<=params.R1 ) = 1.0;
    
    % for boundary conditions
    params.mask ( (R>=params.R2)&(R<=params.R3) ) = 1.0;
end