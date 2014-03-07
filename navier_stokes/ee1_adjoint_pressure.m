function [u_new, uk_tilde, pk] = EE1_adjoint_pressure(u,uk,pk,dt)
    global params
    pk = uk*0; % dummy argument    
    
    nlk = rhs_up (uk,u,'no','yes');
    uk_tilde = uk + dt*nlk;
    uk_tilde = dealias_2d( uk_tilde );    
      
    divubar = divergence_2d( uk_tilde );
      
      for i = 1: 300
          divutilde = divubar + divergence_2d( pk );          
          uk_star = gradient_2d( divutilde );
          testmode = cofitxy(divergence_2d(uk_star));
          
          s= - 0.1*( sum(sum( testmode.*cofitxy(divutilde)))  )/ (sum(sum(testmode.*testmode)) );
          pk = pk + s*uk_star;
%           uk_tilde = uk + dt*nlk - pk;
          
divut(i) = max(max(abs(cofitxy(divutilde))));
          fprintf('field divergence: %12.5e\n',divut(i))
          
          figure(1)
          pcolor(cofitxy(divutilde))
          shading flat
          colorbar
          drawnow
      end
      
      figure(20)
      plot(divut)
        
    u_new = cofitxy_2d(uk_tilde);    
    error('K')
end

% for  n = 1:25
% divUtilde = divUbar + Div(params,p_x,p_y) ;
% disp(max(max( divUtilde )))
% [u_star, v_star] = Grad(params, divUtilde) ;
% testMode = Div(params,u_star,v_star) ;
% 
% s= - ( sum(sum(testMode.*divUtilde))  )/ (sum(sum(testMode.*testMode )) );
% p_x = p_x + s* u_star ;
% p_y = p_y + s* v_star ;
% 
% end
