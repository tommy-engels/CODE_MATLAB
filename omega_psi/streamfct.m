function stream=streamfct(vork)
    global params
    % rotational part of streamfunction
    psi = -cofitxy(poisson(vork));
    
    % potential part depends on mean flow forcing
    str = zeros(params.nx,params.ny);
    switch params.meanflow
        case 'impulsively_x'
            u_inf = 1.0;
            for iy=1:params.ny
                str(:,iy) = u_inf*params.y(iy);
            end
        otherwise
            error('params.meanflow undefined')
    end
   
    % complete streamfunction is both contributions
    stream = psi + str;
end    