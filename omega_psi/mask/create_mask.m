function create_mask
    global params

    %% allocate
    params.mask = zeros(params.nx,params.ny);
    params.us   = zeros(params.nx,params.ny,2);
    params.u_BC = zeros(params.nx,params.ny,2);
    
    
    switch params.imask
        case 'guermond'
            a=-1;
            b=+1;
            [dummy,i]=min(abs(params.x-a));
            [dummy,j]=min(abs(params.x-b));
            
            params.mask = ones(params.nx,params.ny);
            params.mask ( i+1:j-1,i+1:j-1 ) = 0;
            
        case 'empty'
            params.mask = zeros( params.nx,params.ny );
            
        otherwise
            error('params.imask not set')
    end
    
    
    
    if (sum(sum(params.mask))==0)
        warning('Mask is empty..')
    end
end