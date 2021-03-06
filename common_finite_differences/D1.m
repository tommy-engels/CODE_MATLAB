function D=D1(N,h)
    % Returns a derivative matrix 
    % non-periodic
    % Second Order
    mitte=zeros(1,N);
    oben =ones(1,N-1);
    unten=(-1).*ones(1,N-1);

    D=diag(mitte)+diag(oben,1)+diag(unten,-1);
    
    D(1,1)=-3;
    D(1,2)=4;
    D(1,3)=-1;

    D(N,N-2)=1;
    D(N,N-1)=-4;
    D(N,N)=3;
 
    D = sparse(D) /(2*h);  
end 
