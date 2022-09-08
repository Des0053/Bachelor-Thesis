s=2;
a=1;
N=fix((tfin-t0)/dt+1);
M=fix((g_fin-g_in)/dx)+(N-1)*s+1;
z=zeros(N,M);
y=zeros(N,M);
j=1.;
for g0=g_in:dg:g_fin
g(2,j)=g0;
f(2,j)= dt*N0*g0^(-p);
j=j+1;
end
i=2;
for t=t0+dt:dt:tfin
j=1;
while j<M+1 && g(i,j)> 0.0
      g(i+1,j)=g(i,j)+dt*a*g(i,j)/(2*t);
    if g(i,j)>=g_min && g(i,j)<=g_max
     f(i+1,j)=f(i,j)+dt*(N0*g(i+1,j)^(-p)-a*f(i,j)/(2*t));
    else
     f(i+1,j)=f(i,j)-dt*a*f(i,j)/(2*t);
    end
    j=j+1;  
end
j=j-1;
if min(g(i+1,1:j))>g_min
    dx=(min(g(i+1,1:j))-g_min)/(s-1);
for j=j:j+s-1
            g(i+1,j+1)=g_min+dz;
            f(i+1,j+1)=dt*N0*g(i+1,j+1)^(-p);
end
end
i=i+1;
end