function[xvel,yvel,Press,dummyu]=init(N,y)
%Initializing U,V, P
kl=0;
xd=zeros(N);
xvel1=zeros(N);
yvel=zeros(N);
Press=zeros(N);
for j=((N+1)/2)+1:N
    xd(1,j)=1-(y(j,1)-1)^2;
end
for j=1:(N+1)/2
    xd(1,j)= xd(1,(N)-kl);
    kl=kl+1;
end
xd(1,(N+1)/2)=1;
for j=(N+1)/2:N
    xvel1(1,j)=xd(1,2*j-N);
end
xvel=xvel1';
dummyu=xvel;
end
