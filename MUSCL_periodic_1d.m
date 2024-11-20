
N=255*9;
x = [-1:2/(255*9):1];
dx=x(2)-x(1);
%u0_int=cos(x*pi)/(pi*dx)+0.5*x/dx;
%u0=u0_int(2:end)-u0_int(1:end-1);
u0=zeros(1,255*9);
u0(N/4:N/2)=1;
%u0(1:67)=1;
%u0(90:130)=-0.25;
%u0(160:200)=0.5;
x0=[-1+1/(255*9):2/(255*9):1-1/(255*9)];
FV=zeros(100,(255*9));
for i=1:100
    t=(i-1)*0.01;
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
%%
x=[-1+1/255:2/255:1-1/255];
t=[0.01:0.01:1];
fv=surf(x,t,FV(1:100,5:9:end))
fv.EdgeColor = 'none';
%%
usol=FV(1:100,5:9:end)';
save("movingsquare_exact.mat",'t','x','usol')
%%
aa=2*randi(150,1,50);
bb=2*randi(150,1,50);
cc=2*randi(122,1,50);
dd=2*randi(122,1,50);
difference=0;
for i =1:50
    a=aa(i);b=bb(i);c=cc(i);d=dd(i);
    sum_u_dx=(dx/3)*((FV(b,c)-FV(a,c))+4*(sum(FV(b,c+1:2:d-1)-FV(a,c+1:2:d-1)))+2*(sum(FV(b,c+2:2:d-2)-FV(a,c+2:2:d-2)))+(FV(b,d)-FV(a,d)));
    f=FV.*FV/2;
    sum_f_dt=(dx/3)*((f(a,d)-f(a,c))+4*(sum(f(a+1:2:b-1,d)-f(a+1:2:b-1,c)))+2*(sum(f(a+2:2:b-2,d)-f(a+2:2:b-2,c)))+(f(b,d)-f(b,c)));
    difference=difference+sum_u_dx+sum_f_dt;
end
difference/50
%%
diff=0;
aaa=2*randi(150,1,100);
ccc=2*randi(122,1,100);
aaa=sort(aaa);
ccc=sort(ccc);
aa=aaa(1:50);
bb=aaa(51:100);
cc=ccc(1:50);
dd=ccc(51:100);
for i =1:50
    a=aa(i);b=bb(i);c=cc(i);d=dd(i);
    sum_u_dx=(dx/24)*(9*(FV(b,c)-FV(a,c)+FV(b,d)-FV(a,d))+28*(FV(b,c+1)-FV(a,c+1)+FV(b,d-1)-FV(a,d-1))+23*(FV(b,c+2)-FV(a,c+2)+FV(b,d-2)-FV(a,d-2))+24*(sum(FV(b,c+3:d-3)-FV(a,c+3:d-3))));
    f=FV.*FV/2;
    sum_f_dt=(dx/24)*(9*(f(a,d)-f(a,c)+f(b,d)-f(b,c))+28*(f(a+1,d)-f(a+1,c)+f(b-1,d)-f(b-1,c))+23*(f(a+2,d)-f(a+2,c)+f(b-2,d)-f(b-2,c))+24*(sum(f(a+3:b-3,d)-f(a+3:b-3,c))));
    diff=diff+sum_u_dx+sum_f_dt;
end
diff/50

%%
N=256;
x = [-1:2/(7*255):1];
dx=x(2)-x(1);
u0_int=cos(x*pi)/(pi*dx)+0.5*x/dx;
u0=u0_int(2:end)-u0_int(1:end-1);
x0=[-1+1/(255*7):2/(255*7):1-1/(7*255)];
FV=zeros(2101,7*255);
for i=1:2101
    t=(i-1)*dx;
   [t,u]=ode45(@ddt_fv,[0,dx],u0);
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

%%
insum_u_dx=zeros(2101,255);
for i=1:255
    insum_u_dx(:,i)=(dx/24)*(9*(FV(:,1)+FV(:,7*i))+28*(FV(:,2)+FV(:,7*i-1))+23*(FV(:,3)+FV(:,7*i-2))+24*sum(FV(:,4:7*i-3),2));
end
%coarse=zeros(301,255);
coarse=insum_u_dx(1:7:2101,:);
x=[-1+1/(255):2/(255):1-1/(255)];
t=[0:dx:200*dx];
fv=surf(x,t,upred)
fv.EdgeColor = 'none';
%%

%cfv=FV(1:7:end,1:7:end);

usol=usol';
x=[-1+1/(255):2/(255):1-1/(255)];
t=[0:dx:200*dx];
fv=surf(x,t,upred)
fv.EdgeColor = 'none';


%%
x=[-1+1/(255):2/(255):1-1/(255)];
t=[0:dx:200*dx];
usol=FV';
save("step_255_201.mat",'t','x','usol')
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

dudtL=(w-[w(end);w(1:end-1)])/dx;
dudtR=([w(2:end);w(1)]-w)/dx;

s=minmod(dudtR,dudtL);

w_R=w+s*dx/2;
w_L=w-s*dx/2;

fR = Godunov(w_R, [w_L(2:end);w_L(1)]);

fL = [fR(end);fR(1:end-1)];

dwdt = ((fL-fR)/dx);
end

function f = Godunov(wL,wR)
fL = wL.^2/2;
fR = wR.^2/2;
f=max(fR,fL);
f(wL<wR) = min(fL(wL<wR),fR(wL<wR));
f((wL<0)&(wR>0))=0;
end

function s = minmod(uL,uR)
s=min(uL,uR);
s((uL<0) & (uR<0))=max(uL((uL<0) & (uR<0)),uR((uL<0) & (uR<0)));
s(uL.*uR<=0)=0;
end