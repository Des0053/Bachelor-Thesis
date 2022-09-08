for t=2:fix(N/10):N
l=1; i=1;
for omega=lowlim:dw1:uplim
    for j=1:length(gamma(t,:))
        if z(t,j)>0
        omega_c(j)=gamma(t,j)^2*q*B0*(tval(t)/t0)^a/
        (mass*c_light)*sin(a_ph);
        x(i,j)=10^(omega)/omega_c(j);
        end
    end
            i=i+1;
end
zmax=log10(100*max(max(x)));
i=1;
for omega=lowlim:dw1:uplim
    for j=1:length(z(t,:))-1
        if z(t,j)>0 && z(t,j+1)>0
         F(i,j)=Ffun(x(i,j),zmax);
        Ptot(i,l)=Ptot(i,l)+F(i,j)*f(t,j)*sqrt(3)*q^3*B0*(tval(t)/t0)^a*
        sin(a_ph)/(2*pi*mass*c_light^2)*log(10)*gamma(t,j)*
        (log10(gamma(t,j+1))-log10(gamma(t,j)));
        end
        w(i,1)=10^omega;
        freq(i,1)=w(i,1)/(2*pi);
        Pv_tot(i,l)=Ptot(i,l)*freq(i,1);
    end
    i=i+1;
end
l=l+1
end