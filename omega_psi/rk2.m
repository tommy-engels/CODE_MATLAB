function [vork_new,dt] = rk2(time,vork)
    global params  
    % compute non-linear terms
    [nlk,dt] = rhs_omega(time,vork,'yes');
    % advance in time
    vis = exp(-params.nu*dt*(params.Kx.^2+params.Ky.^2) );    
    vork_new = vis.*(vork + dt*nlk );         
    
    % end of euler step    
    [nlk2,dt_dummy] = rhs_omega(time+dt,vork_new,'yes');     
    % advance in time
    vork_new = vis.*vork + 0.5*dt*(nlk.*vis + nlk2);  
    % Dealiasing
    vork_new = dealias(vork_new);
end