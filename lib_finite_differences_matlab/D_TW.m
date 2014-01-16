function D = D_TW(N,h)
    %------------------------------------------------------------- 
    % Expicit Derivate Matrix with improved wave numbers of Tam and Webb
    % Periodic Version
    % D = D_TW(N,h)
    % N number of Points
    % h = Delta x   
    %-------------------------------------------------------------     

    alpha2=-1/(h)*0.18941;
    alpha3= 1/(h)*0.02652; 
    alpha1= 1/(2*h)-2*alpha2-3*alpha3;

    alpha1Old= 1/(h)*0.79927; 
 
    D= alpha1*(diag(ones(N-1,1),1)-  diag(ones(N-1,1),-1)  + diag(ones(1,1),-N+1) -  diag(ones(1,1),N-1))    +...  
    alpha2*(diag(ones(N-2,1),2)-  diag(ones(N-2,1),-2) +  diag(ones(2,1),-(N-2))- diag(ones(2,1),N-2))  +...        
    alpha3*(diag(ones(N-3,1),3)-  diag(ones(N-3,1),-3) +  diag(ones(3,1),-(N-3))- diag(ones(3,1),N-3) );
end