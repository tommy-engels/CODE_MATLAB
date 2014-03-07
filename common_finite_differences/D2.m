
function d=D2(N,h)
    % Returns a derivative matrix 
    % non-periodic
    % Second Order
    mitte=(-2)*ones(1,N);
    oben =ones(1,N-1);
    unten=ones(1,N-1);
    
    d=diag(mitte)+diag(oben,1)+diag(unten,-1);
    
    d(1,1)=2;
    d(1,2)=-5;
    d(1,3)=4;
    d(1,4)=-1;
    
    d(N,N-3)=-1;
    d(N,N-2)=4;
    d(N,N-1)=-5;
    d(N,N)=2;
    
    d= sparse(d) / (h^2);
end 