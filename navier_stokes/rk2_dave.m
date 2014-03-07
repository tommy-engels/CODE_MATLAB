function [u_new, uk_new, pk] = RK2_dave(u,uk,pk,dt)
%% Dave's scheme with modified pressure poisson eqn
    global params
    pk = pk*0; % dummy argument
    
    % compute non-linear terms + modified pressure
    nlk = rhs_up (uk, u, 'no');
    divu = cofitxy(divergence_2d(uk));
    qk = fft2( ((1-params.mask)/params.eta).*divu ) + divergence_2d(nlk);
    qk = poisson(qk);
    % add penalization term, possibly with active penalization
    params.us = create_us( u );
    penal(:,:,1) = fft2((params.mask/params.eta).*(u(:,:,1)-params.us(:,:,1)));
    penal(:,:,2) = fft2((params.mask/params.eta).*(u(:,:,2)-params.us(:,:,2)));
    nlk = nlk - gradient_2d(qk) - penal;
    
    % exponential factor:
    vis = cal_vis_2d( dt );
    % euler step:
    uk_new = vis .* (uk + dt*nlk );
    uk_new = dealias_2d( uk_new );
    u_new = cofitxy_2d (uk_new);
        
    % end of euler step , now doing correction step
    nlk2 = rhs_up(uk_new, u_new,  'no');    
    divu = cofitxy(divergence_2d(uk_new));
    qk = fft2( ((1-params.mask)/params.eta).*divu ) + divergence_2d(nlk2);
    qk = poisson(qk);
    % add penalization term, possibly with active penalization
    params.us = create_us( u );
    penal(:,:,1) = fft2((params.mask/params.eta).*(u_new(:,:,1)-params.us(:,:,1)));
    penal(:,:,2) = fft2((params.mask/params.eta).*(u_new(:,:,2)-params.us(:,:,2)));
    nlk2 = nlk2 - gradient_2d(qk) - penal;
    
    % advance in time
    uk_new = vis.*uk + 0.5*dt*(nlk.*vis + nlk2);
    uk_new = dealias_2d( uk_new );
        
    % go to physical space
    u_new = cofitxy_2d(uk_new);
    
    if (max(max(max(abs(u))))>1e2)
        error('RK2_dave diverged..')
    end
end
