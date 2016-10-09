function[rcs,rus,rvs]=RHS(N,JC,xvel,yvel,Press,zx,ey,ex,zy,Re,ep,eps)
U=zeros(N);
V=zeros(N);
g11=zeros(N);
g12=zeros(N);
g22=zeros(N);
Estar1=zeros(N,N,3);
Estar2=zeros(N,N,3);
Estarv1=zeros(N,N,3);
Estarv2=zeros(N,N,3);
dEs1=zeros(N,N,3);
dEs2=zeros(N,N,3);
dEsv1=zeros(N,N,3);
dEsv2=zeros(N,N,3);
rcs=zeros(N);
rus=zeros(N);
rvs=zeros(N);
JChx=zeros(N);
JChy=zeros(N);
zxhx=zeros(N);
zyhx=zeros(N);
exhx=zeros(N);
eyhx=zeros(N);
zxhy=zeros(N);
zyhy=zeros(N);
exhy=zeros(N);
eyhy=zeros(N);
g11hx=zeros(N);
g11hy=zeros(N);
g12hx=zeros(N);
g12hy=zeros(N);
g22hx=zeros(N);
g22hy=zeros(N);
for i=1:N
    for j=1:N
        U(i,j)=(xvel(i,j)*zx(i,j)+yvel(i,j)*zy(i,j));
        V(i,j)=(xvel(i,j)*ex(i,j)+yvel(i,j)*ey(i,j));
    end
end
for i=1:N
    for j=1:N
        Estar1(i,j,1)=U(i,j)/JC(i,j);
        Estar1(i,j,2)=(xvel(i,j)*U(i,j)+Press(i,j)*zx(i,j))/JC(i,j);
        Estar1(i,j,3)=(yvel(i,j)*U(i,j)+Press(i,j)*zy(i,j))/JC(i,j);
        Estar2(i,j,1)=V(i,j)/JC(i,j);
        Estar2(i,j,2)=(xvel(i,j)*V(i,j)+Press(i,j)*ex(i,j))/JC(i,j);
        Estar2(i,j,3)=(yvel(i,j)*V(i,j)+Press(i,j)*ey(i,j))/JC(i,j);
    end
end
for k=1:3
    for i=2:N-1
        for j=2:N-1
            dEs1(i,j,k)=0.5*(Estar1(i,j+1,k)-Estar1(i,j-1,k));
            dEs2(i,j,k)=0.5*(Estar2(i+1,j,k)-Estar2(i-1,j,k));
        end
    end
end
for i=1:N
    for j=1:N
        g11(i,j)=((zx(i,j))^2)+((zy(i,j))^2);
        g12(i,j)=zx(i,j)*ex(i,j)+zy(i,j)*ey(i,j);
        g22(i,j)=((ex(i,j))^2)+((ey(i,j))^2);
    end
end
for i=1:N-1
    for j=1:N
        JChy(i,j)=0.5*(JC(i+1,j)+JC(i,j));
        zxhy(i,j)=0.5*(zx(i+1,j)+zx(i,j));
        zyhy(i,j)=0.5*(zy(i+1,j)+zy(i,j));
        exhy(i,j)=0.5*(ex(i+1,j)+ex(i,j));
        eyhy(i,j)=0.5*(ey(i+1,j)+ey(i,j));
    end
end
for j=1:N
    JChy(N,j)=0.5*(JC(N-1,j)+JC(N,j));
    zxhy(N,j)=0.5*(zx(N-1,j)+zx(N,j));
    zyhy(N,j)=0.5*(zy(N-1,j)+zy(N,j));
    exhy(N,j)=0.5*(ex(N-1,j)+ex(N,j));
    eyhy(N,j)=0.5*(ey(N-1,j)+ey(N,j));
end
for i=1:N
    for j=1:N-1
        JChx(i,j)=0.5*(JC(i,j+1)+JC(i,j));
        zxhx(i,j)=0.5*(zx(i,j+1)+zx(i,j));
        zyhx(i,j)=0.5*(zy(i,j+1)+zy(i,j));
        exhx(i,j)=0.5*(ex(i,j+1)+ex(i,j));
        eyhx(i,j)=0.5*(ey(i,j+1)+ey(i,j));
    end
end
for i=1:N
    JChx(i,N)=0.5*(JC(i,N-1)+JC(i,N));
    zyhx(i,N)=0.5*(zy(i,N-1)+zy(i,N));
    zxhx(i,N)=0.5*(zx(i,N-1)+zx(i,N));
    exhx(i,N)=0.5*(ex(i,N-1)+ex(i,N));
    eyhx(i,N)=0.5*(ey(i,N-1)+ey(i,N));
end
for i=1:N
    for j=1:N
        g11hx(i,j)=((zxhx(i,j))^2)+((zyhx(i,j))^2);
        g12hx(i,j)=zxhx(i,j)*exhx(i,j)+zyhx(i,j)*eyhx(i,j);
        g22hx(i,j)=((exhx(i,j))^2)+((eyhx(i,j))^2);
        g11hy(i,j)=((zxhy(i,j))^2)+((zyhy(i,j))^2);
        g12hy(i,j)=zxhy(i,j)*exhy(i,j)+zyhy(i,j)*eyhy(i,j);
        g22hy(i,j)=((exhy(i,j))^2)+((eyhy(i,j))^2);
    end
end
for i=1:N-1
    for j=1:N-1
        Estarv1(i,j,1)=0;
        Estarv1(i,j,2)=(g11hx(i,j)*(xvel(i,j+1)-xvel(i,j))+g12hx(i,j)*(xvel(i+1,j)-xvel(i,j)))/(Re*JChx(i,j));
        Estarv1(i,j,3)=(g11hx(i,j)*(yvel(i,j+1)-yvel(i,j))+g12hx(i,j)*(yvel(i+1,j)-yvel(i,j)))/(Re*JChx(i,j));
    end
end
for i=1:N-1
    for j=1:N-1
        Estarv2(i,j,1)=0;
        Estarv2(i,j,2)=(g12hy(i,j)*(xvel(i,j+1)-xvel(i,j))+g22hy(i,j)*(xvel(i+1,j)-xvel(i,j)))/(Re*JChy(i,j));
        Estarv2(i,j,3)=(g12hy(i,j)*(yvel(i,j+1)-yvel(i,j))+g22hy(i,j)*(yvel(i+1,j)-yvel(i,j)))/(Re*JChy(i,j));
    end
end
for i=2:N-1
    for j=2:N-1
        dEsv1(i,j,1)=(Estarv1(i,j,1)-Estarv1(i,j-1,1));
        dEsv1(i,j,2)=(Estarv1(i,j,2)-Estarv1(i,j-1,2));
        dEsv1(i,j,3)=(Estarv1(i,j,3)-Estarv1(i,j-1,3));
        dEsv2(i,j,1)=(Estarv2(i,j,1)-Estarv2(i-1,j,1));
        dEsv2(i,j,2)=(Estarv2(i,j,2)-Estarv2(i-1,j,2));
        dEsv2(i,j,3)=(Estarv2(i,j,3)-Estarv2(i-1,j,3));
    end
end
for i=2:N-1
    for j=2:N-1
        rcs(i,j)=(dEs1(i,j,1)+dEs2(i,j,1)-dEsv1(i,j,1)-dEsv2(i,j,1))*JC(i,j);
        rus(i,j)=(dEs1(i,j,2)+dEs2(i,j,2)-dEsv1(i,j,2)-dEsv2(i,j,2))*JC(i,j);
        rvs(i,j)=(dEs1(i,j,3)+dEs2(i,j,3)-dEsv1(i,j,3)-dEsv2(i,j,3))*JC(i,j);
    end
end
[rho1,rho2]=spectralradius(N,JC,U,V,g11,g22);
[rcs,rus,rvs]=Diss(N,eps,rho1,rho2,Press,xvel,yvel,rcs,rus,rvs,JC);%Addingdissipation
[rcs,rus,rvs]=smoothing(N,ep,rcs,rus,rvs);%addingresidualsmoothing
end
