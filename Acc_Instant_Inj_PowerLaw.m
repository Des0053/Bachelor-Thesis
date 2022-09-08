tval=t0:dt:tfin;
N=length(tval);
M=fix((g_max-g_min)/dg)+1;
a=1;
c=@(g,t,f) g*a*(0.5/t);
u=@(g,t,f)  N0*g^(-p)*exp(-((t-t0)/0.001)^2)-a*(0.5/t)*f;
f0=@(g) N0*g^(-p);
i=1;
for g0=g_min:dg:g_max
    G= @(t,w) [c(w(1), t, w(2)) ; u(w(1), t, w(2))];
    [ts,ws] = ode45(G,tval,[g0;f0(g0)]);
    gamma(:,i)=ws(:,1); 
    dist_fun(:,i)=ws(:,2);
    i=i+1;
end