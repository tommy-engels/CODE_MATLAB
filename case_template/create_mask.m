function create_mask
    global params

    % allocate
    params.mask = zeros(params.nx, params.ny);
    params.us   = zeros(params.nx,params.ny,2);
end
