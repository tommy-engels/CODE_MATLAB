%% this is a wrapper for the active penalization algorithms
function us = create_us( u )
global params

if strcmp(params.active,'passive')
    % classic penalization
    us = params.us; % params.us contains inhomogeneous BC field
elseif strcmp(params.iactive,'chantalat')
    u_smooth = smooth_extension_chantalat ( u, iActive );
    if strcmp(params.CASE,'lamballais')
        % delete smooth extension in outer cylinder 
        % only for the lamballais case
        uu=params.us(:,:,1);
        uu(params.R>params.R1*1.25) = 0;
        params.us(:,:,1)=uu;
        uu=params.us(:,:,2);
        uu(params.R>params.R1*1.25) = 0;
        params.us(:,:,2)=uu;
    end    
    us  = u_smooth + params.u_BC;
elseif strcmp(params.iactive,'dave')
    
end

end