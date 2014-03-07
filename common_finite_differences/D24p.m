function D=D24p(n,h)
    D=   diag(ones(n,1)*(-5/2)) ...
        +diag(ones(n-1,1)*(4/3),1) ...
        +diag(ones(n-1,1)*(4/3),-1)...
        +diag(ones(n-2,1)*(-1/12),2)...
        +diag(ones(n-2,1)*(-1/12),-2);

    D(1,n)=4/3;
    D(1,n-1) = -1/12;
    D(2,n)=-1/12;

    D(n,1)=4/3;
    D(n,2) = -1/12;
    D(n-1,1)=-1/12;

    D=sparse(D/(h^2));
end
