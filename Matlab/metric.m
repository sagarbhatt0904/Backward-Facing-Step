function[JC,zx,ey,zy,ex]=metric(N,x,y)
xz=zeros(1,N);
ye=zeros(1,N);
G=zeros(N);
JC=zeros(N);
zx=zeros(N);
ey=zeros(N);
zy=zeros(N);
ex=zeros(N);
xz=zeros(N);
ye=zeros(N);
yz=zeros(N);
xe=zeros(N);
%Metriccalculation
for i=2:N-1
    for j=1:N
        xe(i,j)=0.5*(x(i+1,j)-x(i-1,j));
        ye(i,j)=0.5*(y(i+1,j)-y(i-1,j));
    end
end
for j=1:N
    xe(1,j)=(x(2,j)-x(1,j));
    ye(1,j)=(y(2,j)-y(1,j));
end

for j=1:N
    xe(N,j)=(x(N,j)-x(N-1,j));
    ye(N,j)=(y(N,j)-y(N-1,j));
end

for i=1:N
    for j=2:N-1
        yz(i,j)=0.5*(y(i,j+1)-y(i,j-1));
        xz(i,j)=0.5*(x(i,j+1)-x(i,j-1));
    end
end


for i=1:N
    yz(i,1)=(y(i,2)-y(i,1));
    xz(i,1)=(x(i,2)-x(i,1));
end
for i=1:N
    yz(i,N)=(y(i,N)-y(i,N-1));
    xz(i,N)=(x(i,N)-x(i,N-1));
end


for i=1:N
    for j=1:N
        G(i,j)=xz(i,j)*ye(i,j)-xe(i,j)*yz(i,j);
        JC(i,j)=1/G(i,j);
    end
end

for i=1:N
    for j=1:N
        zx(i,j) = ye(i,j)*JC(i,j);
        ex(i,j) = -yz(i,j)*JC(i,j);
        zy(i,j) = -xe(i,j)*JC(i,j);
        ey(i,j) = xz(i,j)*JC(i,j);
    end
end
end