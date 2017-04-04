clear all;
close all;
N=121;                                      %Grid size
st=1;                                       %stretching ratio
Re=400;
eps=0.01;                                   %Epsilon for dissipation                                  
CFL=0.01;
VN=0.01;
ep=0.01;                                    % Residual smoothing factor
L=20;
[x,y,dy]=gridgen_mex(N,st,L);               % Grid generation
[xvel,yvel,Press,dummyu]=init_mex(N,y);     % Initial conditions
[JC,zx,ey,zy,ex]=metric_mex(N,x,y);         %Metric Calculation
u1=zeros(N);
u2=zeros(N);
u3=zeros(N);
u_new1=zeros(N);
u_new2=zeros(N);
u_new3=zeros(N);
u_new=zeros(N);
v1=zeros(N);
v2=zeros(N);
v3=zeros(N);
v_new1=zeros(N);
v_new2=zeros(N);
v_new3=zeros(N);
v_new=zeros(N);
p1=zeros(N);
p2=zeros(N);
p3=zeros(N);
p_new1=zeros(N);
p_new2=zeros(N);
p_new3=zeros(N);
p_new=zeros(N);
ps=1;
dj=1;
while dj>=10^(-6)
    for i=1:N
        for j=1:N
            
            U(i,j) = (xvel(i,j)*zx(i,j) + yvel(i,j)*zy(i,j));
            V(i,j) = (xvel(i,j)*ex(i,j) + yvel(i,j)*ey(i,j));
        end
    end
    for i=1:N
        for j=1:N
            g11(i,j)  = ((zx(i,j))^2) + ((zy(i,j))^2);
            g12(i,j)  = zx(i,j)*ex(i,j) +zy(i,j)*ey(i,j);
            g22(i,j)  = ((ex(i,j))^2) + ((ey(i,j))^2);
        end
    end
    [rho1,rho2]=spectralradius_mex(N,JC,U,V,g11,g22);  % calculating spectral radius
    
    for i=1:N
        for j=1:N
             dtau(i,j)=min(CFL*(1/JC(i,j))/(max(rho1(i,j),rho2(i,j))),(Re*VN/max(g11(i,j),g22(i,j))));
        end
    end
    
    [rcs,rus,rvs]=RHS_mex(N,JC,xvel,yvel,Press,zx,ey,ex,zy,Re,ep,eps);
    for i=2:N-1
        for j=2:N-1
            p_new1=xvel+0.25*(dtau(i,j)*rcs(i,j));
            u_new1=xvel+0.25*(dtau(i,j)*rus(i,j));
            v_new1=xvel+0.25*(dtau(i,j)*rvs(i,j));
        end
    end
    [u_new1,v_new1,p_new1]=BC_mex(N,u_new1,v_new1,p_new1,dummyu);
    
    [RHSc1,RHSu1,RHSv1]=RHS_mex(N,JC,u_new1,v_new1,p_new1,zx,ey,ex,zy,Re,ep,eps);
    for i=2:N-1
        for j=2:N-1
            p_new2(i,j)=Press(i,j)+0.33*(dtau(i,j)*RHSc1(i,j));
            u_new2(i,j)=xvel(i,j)+0.33*(dtau(i,j)*RHSu1(i,j));
            v_new2(i,j)=yvel(i,j)+0.33*(dtau(i,j)*RHSv1(i,j));
        end
    end
    [u_new2,v_new2,p_new2]=BC_mex(N,u_new2,v_new2,p_new2,dummyu);
    
    [RHSc2,RHSu2,RHSv2]=RHS_mex(N,JC,u_new2,v_new2,p_new2,zx,ey,ex,zy,Re,ep,eps);
    for i=2:N-1
        for j=2:N-1
            p_new3(i,j)=Press(i,j)+0.5*(dtau(i,j)*RHSc2(i,j));
            u_new3(i,j)=xvel(i,j)+0.5*(dtau(i,j)*RHSu2(i,j));
            v_new3(i,j)=yvel(i,j)+0.5*(dtau(i,j)*RHSv2(i,j));
        end
    end
    [u_new3,v_new3,p_new3]=BC_mex(N,u_new3,v_new3,p_new3,dummyu);
    
    [RHSc3,RHSu3,RHSv3]=RHS_mex(N,JC,u_new3,v_new3,p_new3,zx,ey,ex,zy,Re,ep,eps);
    for i=2:N-1
        for j=2:N-1
            p_new(i,j)=Press(i,j)+(dtau(i,j)*RHSc3(i,j));
            u_new(i,j)=xvel(i,j)+(dtau(i,j)*RHSu3(i,j));
            v_new(i,j)=yvel(i,j)+(dtau(i,j)*RHSv3(i,j));
        end
    end
    [u_new,v_new,p_new]=BC_mex(N,u_new,v_new,p_new,dummyu);
    
    
    dj= norm(xvel-u_new);
    
%   Updating values    
    xvel=u_new;
    yvel=v_new;
    Press=p_new;
   
if rem(ps,1000)==0
   data='data.mat';
   save(data);
end
    ps=ps+1;
end
