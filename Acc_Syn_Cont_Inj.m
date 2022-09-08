a=1;
s=5;
N=fix((tfin-t0)/dt+1);
M=fix((g_max-g_min)/dg)+(N-1)*s+1;
g=zeros(N,M);
f=zeros(N,M);
j=1.;

for g0=g_min:dg:g_max
   g(2,j)=g0;
   f(2,j)= dt*N0*g0^(-p);
   j=j+1;
end

i=2;
for t=t0+dt:dt:tfin-dt
j=1;
   while j<M+1 && g(i,j)> 0.0
     if g(i,j)>1.0
        g(i+1,j)=g(i,j)+dt*(g(i,j)/(2*t)-L*(t/t0)^(2*a)*g(i,j)*g(i,j));
     else
        g(i+1,j)=g(i,j)+dt*(g(i,j)/(2*t));
     end
     if g(i,j)>=g_min && g(i,j)<=g_max
       f(i+1,j)=f(i,j)+dt*N0*g(i+1,j)^(-p)
       +dt*(-1/(2*(t))*f(i,j)+2.0*(t/t0)^(2*a)*L*f(i,j)*g(i,j));
     else
       f(i+1,j)=f(i,j)+dt*(-f(i,j)/(2*t)+2.0*L*(t/t0)^(2*a)*f(i,j)*g(i,j));
     end
     if g(i+1,j)<1.0
        g(i+1,j)=1.0;
     end
    
    j=j+1;
end

j=j-1;
if max(g(i+1,:))<g_max
    dz=(g_max-max(g(i+1,:)))/(s-1);
    g(i+1,j+1)=max(g(i+1,:))+dz;
    f(i+1,j+1)=dt*N0*g(i+1,j+1)^(-p);
    j=j+1;
    
    if s>2
    for j=j:j+s-2
            g(i+1,j+1)=g(i+1,j)+dz;
            f(i+1,j+1)=dt*N0*g(i+1,j+1)^(-p);
    end
    end
end

if min(z(i+1,1:j))>g_min
    dz=(min(z(i+1,1:j))-x_in)/(s-1);
    z(i+1,j+1)=g_min;
    y(i+1,j+1)=dt*N0*z(i+1,j+1)^(-p);
    j=j+1;
    
    if s>2
    for j=j:j+s-2
            z(i+1,j+1)=z(i+1,j)+dz;
            y(i+1,j+1)=dt*N0*z(i+1,j+1)^(-p);
    end
    end
end

i=i+1
end