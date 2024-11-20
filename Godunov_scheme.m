
N=500;
x = [0:N-1]/(N);
dx=x(2)-x(1);
u0_int=-cos(x*2*pi)/(2*pi*dx)+0.5*x/dx;
u0=u0_int(2:end)-u0_int(1:end-1);
for t=0:0.01:3
   [t,u]=ode45(@ddt_fv,[0,0.01],u0);
   plot([x(1:end-1);x(2:end)],[u0;u0],"k");
   drawnow;
   pause(0.001);
   u0=u(end,:);
end


%x=[1:100]/100*2*pi;
%u0=sin(x);
%[t,u] = ode45(@ddt_fv,[0,1],u0);
%plot(u')

function dwdt = ddt_fv(~,w)
n=length(w);
dx = 2*pi/n;
fR = Godunov(w, [w(2:end);w(1)]);
fL = [fR(end);fR(1:end-1)];
dwdt = (fL-fR)/dx;
end

function f = Godunov(wL,wR)
fL = wL.^2/2;
fR = wR.^2/2;
f=max(fR,fL);
f(wL<wR) = min(fL(wL<wR),fR(wL<wR));
f((wL<0)&(wR>0))=0;
end