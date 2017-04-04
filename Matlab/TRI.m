% Tridiagonal Matrix Solver
function TRI(ibeg,iend,aa,bb,cc,dd)
for i=ibeg:iend
    r=aa(i)/bb(i-1);
    bb(i)=bb(i)-r*cc(i-1);
    dd(i)=dd(i)-r*dd(i-1);  
end
dd(iend)=dd(iend)/bb(iend);
for i=iend-1:-1:ibeg
    dd(i)=(dd(i)-cc(i)*dd(i+1)/bb(i));
end
return
end
