function[rcs_r,rus_r,rvs_r]=smoothing(N,ep,rcs,rus,rvs)
%     Residual Smoothing
    Dummyc=zeros(N);
    Dummyu=zeros(N);
    Dummyv=zeros(N);
    a=zeros(1,N);
    b=zeros(1,N);
    c=zeros(1,N);
    dc=zeros(1,N);
    du=zeros(1,N);
    dv=zeros(1,N);
    rcs_r=zeros(N);
    rus_r=zeros(N);
    rvs_r=zeros(N);
    a(1)=0;
    b(1)=1+2*ep;
    c(1)=-ep;
    a(N)=-ep;
    b(N)=1+2*ep;
    c(N)=0;
    for i=1:N
        dc(1)=rcs(i,1);
        du(1)=rus(i,1);
        dv(1)=rvs(i,1);
        for j=2:N-1
            a(j)=-ep;
            b(j)=1+2*ep;
            c(j)=-ep;
            dc(j)=rcs(i,j);
            du(j)=rus(i,j);
            dv(j)=rvs(i,j);
        end
        dc(N)=rcs(i,N);
        du(N)=rus(i,N);
        dv(N)=rvs(i,N);
        TRI(2,N,a,b,c,dc);
        TRI(2,N,a,b,c,du);
        TRI(2,N,a,b,c,dv);
        Dummyc(i,:)=dc;
        Dummyu(i,:)=du;
        Dummyv(i,:)=dv;
    end
    for i=1:N
        dc(1)=Dummyc(i,1);
        du(1)=Dummyu(i,1);
        dv(1)=Dummyv(i,1);
        for j=2:N-1
            a(j)=-ep;
            b(j)=1+2*ep;
            c(j)=-ep;
            dc(j)=Dummyc(i,j);
            du(j)=Dummyu(i,j);
            dv(j)=Dummyv(i,j);
        end
        dc(N)=Dummyc(i,N);
        du(N)=Dummyc(i,N);
        dv(N)=Dummyc(i,N);
        TRI(2,N,a,b,c,dc);
        TRI(2,N,a,b,c,du);
        TRI(2,N,a,b,c,dv);
        rcs_r(i,:)=dc;
        rus_r(i,:)=du;
        rvs_r(i,:)=dv;
    end
end