function uk = mean_flow_forcing( uk )
%% sets zero mode of a velocity field to a given value depending on the case
    global params
    factor = params.nx*params.ny;
    
    switch params.meanflow
        case 'impulsively_x'
            uk (1,1,1) = -1.0 * factor;
    end
       
    
end