function [u0, uk] = inicond()
    global params   
            
    u0(:,:,1) = params.u_ex(:,:,1);
    u0(:,:,2) = params.u_ex(:,:,2);
    
    % make it incompressible:
    uk = coftxy_2d(u0);
    uk = make_incompressible_2d (uk);
    u0 = cofitxy_2d(uk);
    
    
    fprintf('initial field divergence: %12.5e\n',max(max(abs(cofitxy(divergence_2d(uk))))));
end
