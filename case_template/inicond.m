function [u0, uk] = inicond()
    global params   
    
    uk = coftxy_2d(u0);
    fprintf('initial field divergence: %12.5e\n',max(max(abs(cofitxy(divergence_2d(uk))))));
end
