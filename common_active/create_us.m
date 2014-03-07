%% this is a wrapper for the active penalization algorithms
function us = create_us( u )
    global params

    %% classic penalization
    if strcmp(params.active,'passive')
        %% classic penalization
        us = params.u_BC; % params.us contains inhomogeneous BC field

    
    elseif strcmp(params.active,'chantalat')
        %% active penalization: chantalat's prolongation
        u_smooth = smooth_extension_chantalat ( u );
        if strcmp(params.CASE,'lamballais')
            % delete smooth extension in outer cylinder 
            % only for the lamballais case
            uu = u_smooth(:,:,1);
            uu(params.R>params.R1*1.25) = 0;
            u_smooth(:,:,1)=uu;
            uu = u_smooth(:,:,2);
            uu(params.R>params.R1*1.25) = 0;
            u_smooth(:,:,2)=uu;
        end    
        us  = u_smooth + params.u_BC;

    elseif strcmp(params.active,'dave')   
        %% active penalization: dave's prolongation
        u_smooth = smooth_extension_dave ( u );
        if strcmp(params.CASE,'lamballais')
            % delete smooth extension in outer cylinder 
            % only for the lamballais case
            uu = u_smooth(:,:,1);
            uu(params.R>params.R1*1.25) = 0;
            u_smooth(:,:,1)=uu;
            uu = u_smooth(:,:,2);
            uu(params.R>params.R1*1.25) = 0;
            u_smooth(:,:,2)=uu;
        end    
        us  = u_smooth + params.u_BC;
        
    else
        error ('params.active is wrong')
    end

end