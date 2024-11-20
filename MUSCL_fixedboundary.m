
N=256;
x_rho =[-1:2/255:1];
x_u =[-1:2/255:1];
x =[-1:2/255:1];
dx=x_rho(2)-x_rho(1);
%u0_int=cos(x*pi)./(pi*dx);
%u0=u0_int(2:end)-u0_int(1:end-1)+0.5;
u0=ones(1,255)*3;
u0(1:63)=-3;
FV=zeros(80,255);

for i=1:80
    t=(i-1)*dx;
   [t,u]=ode45(@ddt_fv,[0,0.01],u0);
   FV(i,:)=u(end,:);
   w=u0';
   dx=2/length(w);
   dutL=(w-[w(end);w(1:end-1)])/dx;
   dutR=([w(2:end);w(1)]-w)/dx;
   s=minmod(dutR,dutL);
   plot([x(1:end-1);x(2:end)],[u0-s'.*dx/2;u0+s'.*dx/2],"k");
   drawnow;
   pause(0.001);
   u0=u(end,:);
end
x=[-1+1/255:2/255:1-1/255];
t=[0:0.01:0.80];
u0=zeros(1,255);
u0(1:63)=1;
FV=cat(1,u0,FV)
fv=surf(x,t,FV)
fv.EdgeColor = 'none';
%%
usol=FV';
save("movingcompound_muscl.mat",'t','x','usol')
%x=[1:100]/100*2*pi;
%u0=sin(x);
%[t,u] = ode45(@ddt_fv,[0,1],u0);
%plot(u')

%x=[1:100]/100*2*pi;
%u0=sin(x);
%[t,u] = ode45(@ddt_fv,[0,1],u0);
%plot(u')

function dwdt = ddt_fv(~,w)
n=length(w);
dx = 2/n;
dudt=(w(2:end)-w(1:end-1))/dx;
dudtL=[(w(1)+3)*255/2;dudt];
dudtR=[dudt;(-3-w(end))*255/2];
s=minmod(dudtR,dudtL);
w_R=w+s*dx/2;
w_L=w-s*dx/2;
f = Lax_f(w_R(1:end-1), w_L(2:end));
%fL = Godunov([0.5;w_R(1:end-1)], w_L);
dwdt = ([1;f]-[f;0])/dx;
end

function f = Lax_f(wL,wR)
fL = 4*wL.^2./(4*wL.^2+(1-wL).^2);
fR = 4*wR.^2./(4*wR.^2+(1-wR).^2);
f=0.5*(fL+fR)+(1/2.55)*(wL-wR);
%f(wL<wR) = min(fL(wL<wR),fR(wL<wR));
%f((wL<0)&(wR>0))=0;
end

function s = minmod(uL,uR)
s=min(uL,uR);
s((uL<0) & (uR<0))=max(uL((uL<0) & (uR<0)),uR((uL<0) & (uR<0)));
s(uL.*uR<=0)=0;
end