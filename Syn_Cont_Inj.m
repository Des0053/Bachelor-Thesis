s=2;
a=1;
N=fix((tfin-t0)/dt+1);
M=fix((g_max-g_min)/dg)+(N-1)*(s-1)+1;
g=zeros(N,M);
f=zeros(N,M);
j=1.;
for g0=g_min:dx:g_max
g(2,j)=g0;
f(2,j)= dt*N0*g0^(-p);
j=j+1;
end
i=2;
for t=t0+dt:dt:tfin
j=1;
while j<M+1 && g(i,j)> 0.0
    if g(i,j)>1.0
        g(i+1,j)=g(i,j)-dt*L*(t/t0)^(2*a)*g(i,j)*g(i,j);
    else 
        g(i+1,j)=g(i,j);
    end
    if g(i,j)>=g_min && g(i,j)<=g_max
     f(i+1,j)=f(i,j)+dt*N0*g(i+1,j)^(-p)
     +dt*(2.0*L*(t/t0)^(2*k)*f(i,j)*g(i,j));
     else
     f(i+1,j)=f(i,j)+dt*(2.0*L*(t/t0)^(2*k)*f(i,j)*g(i,j));
    end
    j=j+1;  
end
j=j-1;
     dz=(g_max-g(i+1,j))/(s-1);
 for j=j:j+s-2
             g(i+1,j+1)=g(i+1,j)+dz;
             f(i+1,j+1)=dt*N0*g(i+1,j+1)^(-p);
 end
i=i+1;
end