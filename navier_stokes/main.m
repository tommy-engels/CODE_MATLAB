function main
    clear all
    close all
    restoredefaultpath;
    addpath(genpath('../common_spectral/'))
    addpath(genpath('../common_active/'))
    addpath(genpath('../common_misc/'))
    addpath(genpath('../common_finite_differences/'))
    addpath(genpath('./mask/'))
    addpath(genpath('./inicond/'))
    
    [e] = time_step ( 1e-3, 128, @RK2_classic, 0.3 );
end

