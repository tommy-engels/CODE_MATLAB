function u_smooth = smooth_extension_chantalat ( u )
    global params
    u_smooth = zeros(params.nx,params.ny,2);
    
    %% compute BETA field
    % the beta field is the "normal" derivative of the velocity field. so its
    % actually d/dn of each component of u. The n-direction is indicated by the
    % normal vectors, defined as the gradient of phi.  
    beta = compute_beta ( u );


    %% prolongate solving advection eqn
    % here we solve an advection diffusion equation numerically, in order to
    % prolongate the beta field from the outside of the obstacle to its
    % inside. The advection-diffusion equation uses the (negative) normal
    % vector field as velocity. The diffusion helps regularizing. Note: the
    % velocity field -n is *not* at all divergence-free. The could be reflected
    % in the adv-diff eqn (see wikipedia)
    CFL = 0.98;
    dx = params.dx;
    umax = 1.0;
    dt = CFL * dx / umax;
    Tend = 0.10;
    nt = round(Tend/dt);

    for i=1:nt
        switch active_chantalat_stepper
            case 'upwind'
                beta(:,:,1) = EE1(beta(:,:,1),dt,@RHS_advection_upwind);
                beta(:,:,2) = EE1(beta(:,:,2),dt,@RHS_advection_upwind);
            case 'central'
                beta(:,:,1) = RK4(beta(:,:,1),dt,@RHS_advection_central);
                beta(:,:,2) = RK4(beta(:,:,2),dt,@RHS_advection_central);
        end
    end
    
    
    if (isnan(beta))
        error('NaN in beta...');
    end
    
    %% construct actual us field prolongation    
    u_smooth(:,:,1) = params.mask .* ( params.phi.*beta(:,:,1) );
    u_smooth(:,:,2) = params.mask .* ( params.phi.*beta(:,:,2) ); 
end







function rhs = RHS_advection_upwind ( field, dt )
% right hand side for the transport eqn. the velocity field is given by the
% normal vectors (times -1) n_x and n_y
    global params
    dx = params.dx;
    dy = params.dy;
    rhs = zeros(params.nx,params.ny);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % vectorized upwind scheme
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    ix = 2:params.nx-1;
    iy = 2:params.ny-1;
    
    aplus = zeros(params.nx,params.ny);
    amins = zeros(params.nx,params.ny);
    uplus = zeros(params.nx,params.ny);
    umins = zeros(params.nx,params.ny);
    
    aplus(ix,iy) = max( -params.n_x(ix,iy), 0 );
    amins(ix,iy) = min( -params.n_x(ix,iy), 0 );
    umins(ix,iy) = (field(ix,iy)-field(ix-1,iy)) / dx;
    uplus(ix,iy) = (field(ix+1,iy)-field(ix,iy)) / dx;
    
    rhs(ix,iy) = rhs(ix,iy) + params.n_x(ix,iy).*(aplus(ix,iy).*umins(ix,iy) - amins(ix,iy).*uplus(ix,iy));
    
    % second component    
    aplus(ix,iy) = max( -params.n_y(ix,iy),0 );
    amins(ix,iy) = min( -params.n_y(ix,iy),0 );
    umins(ix,iy) = (field(ix,iy)-field(ix,iy-1)) / dy;
    uplus(ix,iy) = (field(ix,iy+1)-field(ix,iy)) / dy;
     
    rhs(ix,iy) = rhs(ix,iy) + params.n_y(ix,iy).*(aplus(ix,iy).*umins(ix,iy) - amins(ix,iy).*uplus(ix,iy));
    
%     laplace = Laplacian(field);
%     lambda = 0.5*0.5*dt;
    
    rhs = params.mask2.*(rhs );%+ lambda*laplace);
end

function rhs=RHS_advection_central ( field, dt )
% right hand side for the transport eqn. the velocity field is given by the
% normal vectors (times -1) n_x and n_y
    global params
    lambda = 0.5*0.5*dt; % old value: 0.5*0.5*dt
    laplace = Laplacian(field);    


    % tried to include non-solenoidal velocity field
%     grad_x = (params.D1* (field.*params.n_x));
%     grad_y = ((field.*params.n_y)*params.D1');  
%     rhs = params.mask2.*( grad_x + grad_y + lambda*laplace ); 

    % ignore the fact that the velocity field is not divergence-free
    [grad_x, grad_y] = MyGradient( field );
    rhs = params.mask.*( params.n_x.*grad_x + params.n_y.*grad_y + lambda*laplace ); 
end

function u1 = EE1(u0,dt,rhs)
    u1 = u0 + dt * rhs(u0,dt) ;
end
 
function u1 = RK4(u0,dt,rhs)
    k1=rhs(u0, dt);
    k2=rhs(u0+dt/2* k1,dt);
    k3=rhs(u0+dt/2* k2,dt);
    k4=rhs(u0+dt* k3,dt);
    u1 = u0+dt/6 * (k1 + 2*k2 +2*k3+k4) ;
end






    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % version upwind scheme    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     for ix=2:params.nx-1
%         for iy=2:params.ny-1
%             if (params.mask2(ix,iy))
%                % this is the x-gradient
%                a = -params.n_x(ix,iy);
%                aplus = max(a,0);
%                amins = min(a,0);               
%                umins = (field(ix,iy)-field(ix-1,iy)) / dx;
%                uplus = (field(ix+1,iy)-field(ix,iy)) / dx;
%                rhs(ix,iy) = rhs(ix,iy) + sign*params.n_x(ix,iy)*(aplus*umins - amins*uplus);
% 
%                % now y-gradient
%                a = -params.n_y(ix,iy);
%                aplus = max(a,0);
%                amins = min(a,0);               
%                umins = (field(ix,iy)-field(ix,iy-1)) / dy;
%                uplus = (field(ix,iy+1)-field(ix,iy)) / dy;
%                rhs(ix,iy) = rhs(ix,iy) + sign*params.n_y(ix,iy)*(aplus*umins - amins*uplus);               
%             end
%         end
%     end
%     rhs = rhs.*params.mask2;



%     % euler step
%     rhs1_old = RHS_advection_central( beta(:,:,1),dt );
%     rhs2_old = RHS_advection_central( beta(:,:,2),dt );
%     beta(:,:,1) = beta(:,:,1) + dt*rhs1_old;
%     beta(:,:,2) = beta(:,:,2) + dt*rhs2_old;    
%     % AB2 time loop
%     for i=1:nt-1
%         rhs1_new = RHS_advection_central( beta(:,:,1),dt );
%         rhs2_new = RHS_advection_central( beta(:,:,2),dt ); 
%         beta(:,:,1) = beta(:,:,1) + 1.5*dt*rhs1_new - 0.5*dt*rhs1_old;
%         beta(:,:,2) = beta(:,:,2) + 1.5*dt*rhs2_new - 0.5*dt*rhs2_old;
%         rhs1_old=rhs1_new;
%         rhs2_old=rhs2_new;
%     end  
