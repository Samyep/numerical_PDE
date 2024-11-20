
%x = [-1:2/510:1];
%y = [0:2/510:2];

x0=[-1+1/510:2/510:1-1/510];
y0=[1-1/510:-2/510:-1+1/510];
dx=x0(2)-x0(1);
dy=y0(2)-y0(1);
dt=0.4*dx;
u0=zeros(510,510);
u0(256:378,122:255)=1;
usol=zeros(200,510,510);

for i=1:200
    t=(i-1)*dt;
    %[t,u]=ode45(@ddt_fv,[0,dt],u0);
    usol(i,:,:)=u0;
    u1=u0+dt*ddt_fv(dx,u0);
    u2=u1+dt*ddt_fv(dx,u0);
    u0=0.5*(u0+u2);
    fv=surf(x0,y0,u0);
    fv.EdgeColor = 'none';
    drawnow;
    pause(0.01);
   
end
t=[0:dt:199*dt];
save("burger_2d.mat",'t','x0','y0','usol')

%%
x=[-1+1/1020:2/1020:1-1/1020];
t=[0:0.01:3];
fv=surf(x,t,FV)
fv.EdgeColor = 'none';
%%
usol=u0;
save("2d_burgers_510_test.mat",'x0','y0','usol')
%%


function dwdt = ddt_fv(~,w)
rowlength=length(w(1,:));
colength=length(w(:,1));
dx = 2/rowlength;
dy = 2/colength;
XdudtL=(w-cat(2,w(:,end),w(:,1:end-1)))/dx;
XdudtR=(cat(2,w(:,2:end),w(:,1))-w)/dx;

YdudtD=(w-[w(2:end,:);w(1,:)])/dy;
YdudtU=([w(end,:);w(1:end-1,:)]-w)/dy;
Xs=minmod(XdudtR,XdudtL);
Ys=minmod(YdudtD,YdudtU);

w_R=w+Xs*dx/2;
w_L=w-Xs*dx/2;
w_U=w+Ys*dy/2;
w_D=w-Ys*dy/2;
fR = Lax_f( w_R,cat(2,w_L(:,2:end),w_L(:,1)));
fU = Lax_f( w_U,[w_D(end,:);w_D(1:end-1,:)]);
fL = Lax_f2(w_L, cat(2,w_R(:,end),w_R(:,1:end-1)));%cat(2,fR(:,end),fR(:,1:end-1));
fD = Lax_f2(w_D, [w_U(2:end,:);w_U(1,:)]); %[fU(2:end,:);fU(1,:)];

dwdt = -(fL+fR+fD+fU)/(dy*dx);
%disp(dwdt(1:52,1:52))
end


function f = Lax_f(w_in,w_out)
dx=2/510;
f_in = (w_in.^2)/2;
f_out = (w_out.^2)/2;
f=0.5*(f_in*dx+f_out*dx)+(1/2)*(w_in*dx-w_out*dx);
end

function f = Lax_f2(w_in,w_out)
dx=2/510;
f_in = (w_in.^2)/2;
f_out = (w_out.^2)/2;
f=-0.5*(f_in*dx+f_out*dx)+(1/2)*(w_in*dx-w_out*dx);
end

function s = minmod(uL,uR)
s=min(uL,uR);
s((uL<0) & (uR<0))=max(uL((uL<0) & (uR<0)),uR((uL<0) & (uR<0)));
s(uL.*uR<= 0)= 0;
end
