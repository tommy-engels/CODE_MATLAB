function create_mask
    global params

    % allocate
    params.mask = zeros(params.nx, params.ny);
    params.us   = zeros(params.nx,params.ny,2);
    
    % set open cavity mask
    params.mask ( (params.X<=0.0)|(params.X>=2.0) ) = 1; 
end
