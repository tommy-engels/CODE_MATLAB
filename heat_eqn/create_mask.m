function create_mask
    global params
    params.mask = ones(params.nx,params.ny);
    params.us   = zeros(params.nx,params.ny);
    
    [dummy,i]=min(abs(params.x-0));
    [dummy,j]=min(abs(params.x-2*pi));
    
    params.mask ( i+1:j-1,i+1:j-1 ) = 0;
end