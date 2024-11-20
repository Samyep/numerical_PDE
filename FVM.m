N=500;
x = [1:N-1]/(N);
dx=x(2)-x(1);
u0_int=-cos(x*2*pi)/(2*pi*dx)+0.5*x/dx;
u0=u0_int(2:end)-u0_int(1:end-1);
for t=0:0.01:0.8
   [t,u]=ode45(@burgers,[0,0.01],u0);
   plot([x(1:end-1);x(2:end)],[u0;u0],"k");
   drawnow;
   pause(0.01);
   u0=u(end,:);
end


function dudt = burgers(~,u)
dx=1./length(u);

u_intL_R = u;
u_intL_L = [u(end);u(1:end-1)];
u_intL_avg = (u_intL_L+u_intL_R)/2;
u_iniL = (u_intL_avg<0).*u_intL_R+(u_intL_avg>=0).*u_intL_L;
f=(u_iniL.^2)/2;
dudt=(f-[f(2:end);f(1)])/dx;
end