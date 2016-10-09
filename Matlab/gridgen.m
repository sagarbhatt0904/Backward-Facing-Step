function[x,y,dy]=gridgen(N,st,L)
% Intializing the grid
x1=zeros(1,N);
y1=zeros(1,N);
x=zeros(N);
y=zeros(N);
dy=2/(N-1);
y1=0:dy:2;
if st==1
    dx=L/(N-1);
    x1=0:dx:L;
else
dx=(L*(1-st)/(1-(st^(N-1))));
x1(1)=0;
for i=2:N-1
    x1(i)=x1(i-1)+(st^(i-2))*dx;
end
x1(N)=L;
end
for i=1:N
for j=1:N
x(i,j)=x1(j);
y(i,j)=y1(i);
end
end
end
