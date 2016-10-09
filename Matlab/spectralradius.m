function[rho1,rho2]=spectralradius(N,JC,U,V,g11,g22)
%Computing spectral radius
rho1=zeros(N);
rho2=zeros(N);
for i=1:N
    for j=1:N
        rho1(i,j)=(1/JC(i,j))*(abs(U(i,j))+sqrt(((U(i,j))^2)+g11(i,j)));
        rho2(i,j)=(1/JC(i,j))*(abs(V(i,j))+sqrt(((V(i,j))^2)+g22(i,j)));
    end
end
end
