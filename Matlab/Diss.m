function[rcs,rus,rvs]=Diss(N,eps,rho1,rho2,Press,xvel,yvel,rcs,rus,rvs,JC)
%computing dissipation
D1_c=zeros(N);
D1_u=zeros(N);
D1_v=zeros(N);
D2_c=zeros(N);
D2_u=zeros(N);
D2_v=zeros(N);
Diss_c=zeros(N);
Diss_u=zeros(N);
Diss_v=zeros(N);

for i=2:N-2
    for j=2:N-2
        D1_c(i,j)=eps*(rho1(i,j)+rho1(i+1,j))*0.5*(Press(i+2,j)-3*Press(i+1,j)+3*Press(i,j)-Press(i-1,j));
        D1_u(i,j)=eps*(rho1(i,j)+rho1(i+1,j))*0.5*(xvel(i+2,j)-3*xvel(i+1,j)+3*xvel(i,j)-xvel(i-1,j));
        D1_v(i,j)=eps*(rho1(i,j)+rho1(i+1,j))*0.5*(yvel(i+2,j)-3*yvel(i+1,j)+3*yvel(i,j)-yvel(i-1,j));
        D2_c(i,j)=eps*(rho2(i,j)+rho2(i,j+1))*0.5*(Press(i,j+2)-3*Press(i,j+1)+3*Press(i,j)-Press(i,j-1));
        D2_u(i,j)=eps*(rho2(i,j)+rho2(i,j+1))*0.5*(xvel(i,j+2)-3*xvel(i,j+1)+3*xvel(i,j)-xvel(i,j-1));
        D2_v(i,j)=eps*(rho2(i,j)+rho2(i,j+1))*0.5*(yvel(i,j+2)-3*yvel(i,j+1)+3*yvel(i,j)-yvel(i,j-1));
        Diss_c(i,j)=((D1_c(i+1,j)+D1_c(i,j))*0.5-0.5*(D1_c(i,j)+D1_c(i-1,j)))+((D2_c(i+1,j)+D2_c(i,j))*0.5-0.5*(D2_c(i,j)+D2_c(i,j-1)));
        Diss_u(i,j)=((D1_u(i+1,j)+D1_u(i,j))*0.5-0.5*(D1_u(i,j)+D1_u(i-1,j)))+((D2_u(i+1,j)+D2_u(i,j))*0.5-0.5*(D2_u(i,j)+D2_u(i,j-1)));
        Diss_v(i,j)=((D1_v(i+1,j)+D1_v(i,j))*0.5-0.5*(D1_v(i,j)+D1_v(i-1,j)))+((D2_v(i+1,j)+D2_v(i,j))*0.5-0.5*(D2_v(i,j)+D2_v(i,j-1)));
        
    end
end
for i=1:N
    Diss_c(i,1)=Diss_c(i,2);
    Diss_u(i,1)=Diss_u(i,2);
    Diss_v(i,1)=Diss_v(i,2);
    Diss_c(i,N-1)=Diss_c(i,N-2);
    Diss_u(i,N-1)=Diss_u(i,N-2);
    Diss_v(i,N-1)=Diss_v(i,N-2);
    Diss_c(i,N)=Diss_c(i,N-1);
    Diss_u(i,N)=Diss_u(i,N-1);
    Diss_v(i,N)=Diss_v(i,N-1);
end
for j=1:N
    Diss_c(1,j)=(Diss_c(2,j));
    Diss_u(1,j)=(Diss_u(2,j));
    Diss_v(1,j)=(Diss_v(2,j));
    Diss_c(N-1,j)=Diss_c(N-2,j);
    Diss_u(N-1,j)=Diss_u(N-2,j);
    Diss_v(N-1,j)=Diss_v(N-2,j);
    Diss_c(N,j)=Diss_c(N-1,j);
    Diss_u(N,j)=Diss_u(N-1,j);
    Diss_v(N,j)=Diss_v(N-1,j);
end
for i=1:N
    for j=1:N
        rcs(i,j)=rcs(i,j)+(Diss_c(i,j)*JC(i,j));
        rus(i,j)=rus(i,j)+(Diss_u(i,j)*JC(i,j));
        rvs(i,j)=rvs(i,j)+(Diss_v(i,j)*JC(i,j));
    end
end
end
