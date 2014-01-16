function e = error_ref(vor, u)
    global params
    
    u_abs    = sqrt(u(:,:,1).^2+u(:,:,2).^2);
    u_abs_ex = sqrt(params.u_ex(:,:,1).^2+params.u_ex(:,:,2).^2);
    fluid = zeros(params.nx,params.ny);
    R = sqrt ( (params.X-0.5*params.Lx).^2 + (params.Y-0.5*params.Ly).^2 );
    fluid ( (R>=params.R1)&(R<=params.R2-params.dx) ) = 1;
    e = norm(reshape(fluid.*(u_abs-u_abs_ex),[],1)) / norm(reshape(fluid.*(u_abs_ex),[],1));
end