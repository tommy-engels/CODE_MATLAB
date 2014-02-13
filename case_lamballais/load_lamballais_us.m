function load_lamballais_us()
    global params
    
    x0 = 0.5*params.Lx;
    y0 = 0.5*params.Ly;
    us = zeros(params.nx,params.ny,2);
    R  = sqrt( (params.X-x0).^2 + (params.Y-y0).^2); 
    
    
    %% load data from file or generate it
    file = ['./case_lamballais/us_lamballais_long_' num2str(params.nx) '.mat'];
    if (exist(file)~=2)
        fprintf('generating us field since data unavailable..\n');
        for ix=1:params.nx
            tic
            for iy=1:params.ny
                % fetch data from lamballais paper
                x = params.X(ix,iy) - 0.5*params.Lx;
                y = params.Y(ix,iy) - 0.5*params.Ly;
                
                [uu,vv,pp,vvor]=lamballais(x,y);
                us(ix,iy,1) = uu;
                us(ix,iy,2) = vv;            
            end
            t=toc;
            fprintf('%s\n',secs2hms( (params.nx-ix)*t ) )
        end
        save(file,'us')
    end
    load (file)
    
    
    params.u_ex = us;
    
    %% kill us field in cylinder
    for ix=1:params.nx
        for iy=1:params.ny
            if ( R(ix,iy) <= params.R1 )
                params.us(ix,iy,:)=0;
            end
        end
    end
    
    fprintf('us:  %s   field is ready...\n', file);
end