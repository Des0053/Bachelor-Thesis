a=1;
tval=t0:dt:tfin;
c=@(g,t,f) g*(0.5/t-L*(t/t0)^(2*a)*g);
f=@(g,t,f)  N0*g^(-p)*exp(-((t-t0)/0.001)^2)-(0.5/t-2.0*L*(t/t0)^(2*a)*g)*f;
f0=@(g) N0*g^(-p);
M=fix((g_max-g_min)/dg)+1;
N=length(tval);
i=1;
for g0=g_min:dg:g_max
    G= @(t,w) [c(w(1), t, w(2)) ; u(w(1), t, w(2))];
    [ts,ws] = ode45(G,tval,[g0;f0(g0)]);
    gamma(:,i)=ws(:,1); 
    dist_fun(:,i)=ws(:,2);
    i=i+1;
end