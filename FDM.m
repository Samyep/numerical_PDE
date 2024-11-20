N=2551;
x = [0:N-1]/N;
u0=sin(x*2*pi);
%u0=zeros(1,N);
%u0(N/2:end)=1;
FD=zeros(81,256);
for i=1:81
    t=(i-1)*0.01;
   [t,u]=ode45(@burgers,[0,0.01],u0);
   FD(i,:)=u(end,1:10:end);
   plot(x,u(end,:),"--o")
   drawnow;
   pause(0.02);
   u0=u(end,:);
end
t=linspace(0,0.8,81);
x=linspace(0,1,256);
ax=surf(x,t,FD)
as.EdgeColor='none';
usol=FD';
save("vanilla_burgers.mat",'t','x','usol')
function dudt = burgers(~,u)
N=length(u);
dx=1.0/N;
uL=[u(end);u(1:end-1)];
uR=[u(2:end);u(1)];
dudxL = (u-uL)/(dx);
dudxR = (uR-u)/(dx);
dudx = (u>0).*dudxL+(u<0).*dudxR;
%dudx=(uR-u)/(dx);
d2udx2 = (uR-2*u+uL)/(dx^2);
a=1;
k=0.01/pi;
dudt=-u.*dudx +k*d2udx2;
end

