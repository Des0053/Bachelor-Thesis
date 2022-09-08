a=1;
tval=t0:dt:tfin;
c=@(x,t,u) -L*(t/t0)^(2*a)*x^2;
N=length(tval);
M=fix((g_max-g_min)/dg)+1;
i=1;
for  g0=g_min:dg:g_max
    f=@(x,t,u)  N0*x0^(-p)*exp(-((t-t0)/0.001)^2)+2.0*L*(t/t0)^(2*a)*g*f;
    f0=@(g) N0*g0^(-p); 
    G= @(t,w) [c(w(1), t, w(2)) ; u(w(1), t, w(2))];
    [ts,ws] = ode45(G,tval,[g0;f0(g0)]);
    gamma(:,i)=ws(:,1); 
    dist_fun(:,i)=ws(:,2);
    i=i+1;
end