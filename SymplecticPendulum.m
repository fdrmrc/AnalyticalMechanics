clear all
close all

#Risoluzione dell'equazione del pendolo semplice (non lineare) mostrando l'effettiva conservazione dell'energia per un RungeKutta semiimplicito di ordine2, mentre Eulero esplicito non la conserva

#Si tratta di un sistema conservativo con energia potenziale U(y)=-cos(y)


#y''=-sin(y)
f=@(y) [y(2);-sin(y(1))];
df=@(y)[0,1;-cos(y(1)),0];
y0=[1;2]; #condizione iniziale


ti=-5;			#Explicit Euler
tf=5;
ts=1000;%time steps
dt=(tf-ti)/ts;
yE=NaN(2,ts+1);
yE(:,1)=y0;
for n=1:ts
  yE(:,n+1)=yE(:,n)+dt*f(yE(:,n));
end

g=9.81;
EEE=yE(2,:).^2 -2*cos(yE(1,:));


				#Runge-Kutta semi-implicito ordine 2

a11=(3+sqrt(3))/6;
a12=0;
a21=-sqrt(3)/3;
a22=(3+sqrt(3))/6;
b1=0.5;
b2=0.5;

y=NaN(2,ts+1);
y(:,1)=y0;

for n=1:ts
  yn=y(:,n);
  F1=@(x) x-yn-dt*a11*f(x);
  J1=@(x) eye(2)-dt*a11*df(x);
  s1=yn;
  res=-J1(s1)\F1(s1);
  tol=dt^4;
  while (norm(res,inf)>tol)
    s1+=res;
    res=-J1(s1)\F1(s1);
  end
  s1+=res;


  F2=@(x) x-yn-a21*dt*f(s1)-a22*dt*f(x);
  J2=@(x) eye(2)-a22*dt*df(x);
  s2=yn;
  res=-J2(s2)\F2(s2);
  while(norm(res,inf)>tol)
    s2+=res;
    res=-J2(s2)\F2(s2);
  end
  s2+=res;

  y(:,n+1)=yn+dt*(b1*f(s1)+b2*f(s2));

endfor

ERK=y(2,:).^2 - 2*cos(y(1,:));

t=linspace(ti,tf,ts+1);

figure 1
plot(t,EEE,'b-*',t,ERK,'r-*')
title('Energy plot')
legend('Energy Expl. Euler','EnergyRK2SI')
xlabel('t')
ylabel('Energy')

figure 2
plot(yE(1,:),yE(2,:),'-*',y(1,:),y(2,:),'-*')
title('Diagrammi di fase')
xlabel('x')
ylabel('xDot')
legend('eulero','RKSI2')
