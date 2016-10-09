function[u_new,v_new,p_new]=BC(N,u_new,v_new,p_new,dummyu)
u_new(:,1)=dummyu(:,1);
u_new(1,:)=0;
u_new(N,:)=0;
u_new(:,N)=u_new(:,N-1);

v_new(:,1)=0;
v_new(1,:)=0;
v_new(N,:)=0;
v_new(:,N)=v_new(:,N-1);

p_new(1,:)=p_new(2,:);
p_new(:,1)=p_new(:,2);
p_new(N,:)=p_new(N-1,:);
p_new(:,N)=p_new(:,N-1);

end
