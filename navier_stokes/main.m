% function main
%     clear all
%     close all
% %     restoredefaultpath;
%     addpath(genpath('../common_spectral/'))
%     addpath(genpath('../common_active/'))
%     addpath(genpath('../common_misc/'))
%     addpath(genpath('../common_finite_differences/'))
%     addpath(genpath('./mask/'))
%     addpath(genpath('./inicond/'))
%     
%     [e] = time_step ( @rk2_implicit );
% end



function main
    clear all
    close all
%     restoredefaultpath;
    addpath(genpath('../common_spectral/'))
    addpath(genpath('../common_active/'))
    addpath(genpath('../common_misc/'))
    addpath(genpath('../common_finite_differences/'))
    addpath(genpath('./mask/'))
    addpath(genpath('./inicond/'))
    
%     load yo.mat
    
%     dt = 2.5*logspace(-4,-3,5);
%     dt = fliplr(dt);
    dt=1e-3;
    
    eps1=logspace(-4,-2,20);
    eps2=eps1;logspace(-3,-2,10);
        
    for i = 1:length(eps1)
        e1(i) = time_step ( eps1(i),dt,@rk2_implicit );
        e2(i) = time_step ( eps1(i),dt,@rk2_implicit_old );        
    end
    for i=1:length(eps2)
        e3(i) = time_step ( eps2(i),eps2(i),@rk2_classic );
    end
    
    figure
    loglog(eps1,e1, eps1,e2, eps2,e3)
    save yo.mat
end

