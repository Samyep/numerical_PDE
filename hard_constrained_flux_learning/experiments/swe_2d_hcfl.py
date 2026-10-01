
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch import nn
import matplotlib.pyplot as plt

torch.set_num_threads(2)

G = 9.81
NREF = 48
NCOARSE = 12
FACTOR = NREF // NCOARSE
DT_SNAPSHOT = 8e-4
LAMBDA = DT_SNAPSHOT / (1.0 / NCOARSE)
NSNAP = 8


# ============================================================
# NumPy reference solver: periodic 2D shallow water
# ============================================================
def np_flux_x(U):
    h = U[..., 0]
    mx = U[..., 1]
    my = U[..., 2]
    u, v = mx / h, my / h
    return np.stack([mx, mx*u + 0.5*G*h*h, mx*v], axis=-1)


def np_flux_y(U):
    h = U[..., 0]
    mx = U[..., 1]
    my = U[..., 2]
    u, v = mx / h, my / h
    return np.stack([my, my*u, my*v + 0.5*G*h*h], axis=-1)


def np_rusanov_x(UL, UR):
    fL, fR = np_flux_x(UL), np_flux_x(UR)
    hL, hR = UL[...,0], UR[...,0]
    uL, uR = UL[...,1]/hL, UR[...,1]/hR
    a = np.maximum(np.abs(uL)+np.sqrt(G*hL), np.abs(uR)+np.sqrt(G*hR))
    return 0.5*(fL+fR)-0.5*a[...,None]*(UR-UL)


def np_rusanov_y(UL, UR):
    fL, fR = np_flux_y(UL), np_flux_y(UR)
    hL, hR = UL[...,0], UR[...,0]
    vL, vR = UL[...,2]/hL, UR[...,2]/hR
    a = np.maximum(np.abs(vL)+np.sqrt(G*hL), np.abs(vR)+np.sqrt(G*hR))
    return 0.5*(fL+fR)-0.5*a[...,None]*(UR-UL)


def np_rhs(U, dx):
    Fx = np_rusanov_x(U, np.roll(U, -1, axis=-2))
    Gy = np_rusanov_y(U, np.roll(U, -1, axis=-3))
    return -(Fx-np.roll(Fx,1,axis=-2))/dx - (Gy-np.roll(Gy,1,axis=-3))/dx


def np_ssprk2(U, dt, dx):
    U1 = U + dt*np_rhs(U,dx)
    return 0.5*U + 0.5*(U1 + dt*np_rhs(U1,dx))


def generate_ic(ntraj, N, seed, ood=False):
    rng = np.random.default_rng(seed)
    x = (np.arange(N)+0.5)/N
    X, Y = np.meshgrid(x,x,indexing="xy")
    U = np.zeros((ntraj,N,N,3),dtype=np.float64)

    for j in range(ntraj):
        typ = rng.choice(["smooth","dam","quadrants"], p=[0.45,0.35,0.20])

        if typ == "smooth":
            h0 = rng.uniform(0.9,1.4)
            ah = rng.uniform(0.08,0.35 if ood else 0.22)
            au = rng.uniform(0.05,0.7 if ood else 0.4)
            av = rng.uniform(0.05,0.7 if ood else 0.4)
            kx, ky = int(rng.integers(1,4)), int(rng.integers(1,4))
            phx, phy = rng.uniform(0,2*np.pi,2)
            h = h0 + ah*np.sin(2*np.pi*kx*X+phx)*np.cos(2*np.pi*ky*Y+phy)
            h = np.maximum(h,0.45)
            u = au*np.sin(2*np.pi*kx*X+rng.uniform(0,2*np.pi))
            v = av*np.cos(2*np.pi*ky*Y+rng.uniform(0,2*np.pi))

        elif typ == "dam":
            hb = rng.uniform(0.7,1.2)
            hi = rng.uniform(1.3,2.4 if ood else 1.9)
            h = np.full((N,N),hb)
            cx,cy = rng.uniform(.25,.75,2)
            wx,wy = rng.uniform(.10,.25,2)
            mask=(np.abs(X-cx)<wx)&(np.abs(Y-cy)<wy)
            h[mask]=hi
            u=np.full((N,N),rng.uniform(-.25,.25))
            v=np.full((N,N),rng.uniform(-.25,.25))

        else:
            vals = rng.uniform(0.65 if not ood else 0.4,
                               1.8 if not ood else 2.4, size=4)
            h=np.empty((N,N))
            h[(X<.5)&(Y<.5)] = vals[0]
            h[(X>=.5)&(Y<.5)] = vals[1]
            h[(X<.5)&(Y>=.5)] = vals[2]
            h[(X>=.5)&(Y>=.5)] = vals[3]
            umax=.6 if not ood else 1.0
            u=np.full((N,N),rng.uniform(-umax,umax))
            v=np.full((N,N),rng.uniform(-umax,umax))

        U[j,...,0]=h
        U[j,...,1]=h*u
        U[j,...,2]=h*v
    return U


def generate_trajectory(ntraj, seed, ood=False):
    U = generate_ic(ntraj,NREF,seed,ood)
    dx=1.0/NREF

    def restrict(V):
        return V.reshape(ntraj,NCOARSE,FACTOR,NCOARSE,FACTOR,3).mean(axis=(2,4)).astype(np.float32)

    snaps=[restrict(U)]
    for _ in range(NSNAP-1):
        rem=DT_SNAPSHOT
        while rem>1e-14:
            h=U[...,0]
            u=U[...,1]/h
            v=U[...,2]/h
            c=np.sqrt(G*h)
            rate=np.max((np.abs(u)+c)/dx + (np.abs(v)+c)/dx)
            dt=min(rem,0.30/rate)
            U=np_ssprk2(U,dt,dx)
            if U[...,0].min()<=0:
                raise RuntimeError("Reference lost positivity")
            rem-=dt
        snaps.append(restrict(U))
    return np.stack(snaps,axis=1)


# ============================================================
# Torch physics
# ============================================================
def primitive(U):
    h=U[...,0].clamp_min(1e-6)
    return torch.stack([h,U[...,1]/h,U[...,2]/h],dim=-1)


def t_flux_x(U):
    h=U[...,0].clamp_min(1e-6); mx=U[...,1]; my=U[...,2]
    u,v=mx/h,my/h
    return torch.stack([mx,mx*u+0.5*G*h*h,mx*v],dim=-1)


def t_flux_y(U):
    h=U[...,0].clamp_min(1e-6); mx=U[...,1]; my=U[...,2]
    u,v=mx/h,my/h
    return torch.stack([my,my*u,my*v+0.5*G*h*h],dim=-1)


def t_rusanov_x(U):
    UR=torch.roll(U,-1,dims=-2)
    fL,fR=t_flux_x(U),t_flux_x(UR)
    hL,hR=U[...,0].clamp_min(1e-6),UR[...,0].clamp_min(1e-6)
    uL,uR=U[...,1]/hL,UR[...,1]/hR
    a=torch.maximum(torch.abs(uL)+torch.sqrt(G*hL),torch.abs(uR)+torch.sqrt(G*hR))
    return .5*(fL+fR)-.5*a[...,None]*(UR-U)


def t_rusanov_y(U):
    UR=torch.roll(U,-1,dims=-3)
    fL,fR=t_flux_y(U),t_flux_y(UR)
    hL,hR=U[...,0].clamp_min(1e-6),UR[...,0].clamp_min(1e-6)
    vL,vR=U[...,2]/hL,UR[...,2]/hR
    a=torch.maximum(torch.abs(vL)+torch.sqrt(G*hL),torch.abs(vR)+torch.sqrt(G*hR))
    return .5*(fL+fR)-.5*a[...,None]*(UR-U)


def entropy_variables(U):
    h=U[...,0].clamp_min(1e-6)
    u,v=U[...,1]/h,U[...,2]/h
    return torch.stack([G*h-.5*(u*u+v*v),u,v],dim=-1)


def psi_x(U):
    h=U[...,0].clamp_min(1e-6)
    u=U[...,1]/h
    return .5*G*h*h*u


def psi_y(U):
    h=U[...,0].clamp_min(1e-6)
    v=U[...,2]/h
    return .5*G*h*h*v


def hard_project(F,U,direction):
    UR=torch.roll(U,-1,dims=-2 if direction=="x" else -3)
    a=entropy_variables(UR)-entropy_variables(U)
    b=(psi_x(UR)-psi_x(U)) if direction=="x" else (psi_y(UR)-psi_y(U))
    r=(a*F).sum(-1)-b
    n2=(a*a).sum(-1)
    alpha=torch.zeros_like(r)
    mask=(r>0)&(n2>1e-14)
    alpha[mask]=r[mask]/n2[mask]
    return F-alpha[...,None]*a


def entropy_residual(F,U,direction):
    UR=torch.roll(U,-1,dims=-2 if direction=="x" else -3)
    a=entropy_variables(UR)-entropy_variables(U)
    b=(psi_x(UR)-psi_x(U)) if direction=="x" else (psi_y(UR)-psi_y(U))
    return (a*F).sum(-1)-b


def fv_step(U,Fx,Gy):
    return U-LAMBDA*(Fx-torch.roll(Fx,1,dims=-2)+Gy-torch.roll(Gy,1,dims=-3))


# ============================================================
# Shared orientation-aware learned flux
# ============================================================
class FaceNet(nn.Module):
    def __init__(self, mean, std, width=64):
        super().__init__()
        self.net=nn.Sequential(
            nn.Linear(27,width),nn.Tanh(),
            nn.Linear(width,width),nn.Tanh(),
            nn.Linear(width,3),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean",torch.tensor(mean,dtype=torch.float32))
        self.register_buffer("std",torch.tensor(std,dtype=torch.float32))
        self.register_buffer("scale",torch.tensor([.5,2.0,1.0],dtype=torch.float32))

    def oriented_patch(self,U,direction):
        P=primitive(U)
        if direction=="x":
            Q=P
        else:
            # transpose y/x so normal direction becomes last spatial axis;
            # swap normal/tangential velocities: [h,v,u].
            Q=P.transpose(-3,-2)[...,[0,2,1]]

        patches=[]
        for dy in [-1,0,1]:
            for dx in [-1,0,1]:
                R=torch.roll(Q,(dy,dx),dims=(-3,-2))
                patches.append((R-self.mean)/self.std)
        feat=torch.cat(patches,dim=-1)
        return feat,Q

    def forward(self,U,direction):
        feat,Q=self.oriented_patch(U,direction)
        corr=torch.tanh(self.net(feat))
        QR=torch.roll(Q,-1,dims=-2)
        jump=torch.sqrt((((QR-Q)/self.std)**2).sum(-1)+1e-12)

        if direction=="x":
            base=t_rusanov_x(U)
            F=base + .22*jump[...,None]*corr*self.scale
            return F
        else:
            # network output is [mass, normal-mom(y), tangential-mom(x)]
            base_oriented=t_rusanov_y(U).transpose(-3,-2)[...,[0,2,1]]
            For=base_oriented + .22*jump[...,None]*corr*self.scale
            # map back and transpose spatial axes
            return For[...,[0,2,1]].transpose(-3,-2)


class Solver(nn.Module):
    def __init__(self,stats,mode="plain",soft_weight=1e-2):
        super().__init__()
        self.mode=mode
        self.soft_weight=soft_weight
        self.face_net=FaceNet(*stats)

    def fluxes(self,U):
        Fx=self.face_net(U,"x")
        Gy=self.face_net(U,"y")
        if self.mode=="hard":
            Fx=hard_project(Fx,U,"x")
            Gy=hard_project(Gy,U,"y")
        return Fx,Gy

    def one_step(self,U):
        Fx=self.face_net(U,"x");Gy=self.face_net(U,"y")
        penalty=torch.tensor(0.,dtype=U.dtype)
        if self.mode=="hard":
            Fx=hard_project(Fx,U,"x");Gy=hard_project(Gy,U,"y")
        elif self.mode=="soft":
            penalty=(torch.relu(entropy_residual(Fx,U,"x")).pow(2).mean()
                     +torch.relu(entropy_residual(Gy,U,"y")).pow(2).mean())
        return fv_step(U,Fx,Gy),penalty


def train_seed(seed,outdir="/mnt/data/swe_2d_hcfl_runs",iterations=650):
    out=Path(outdir);out.mkdir(exist_ok=True)
    train_np=generate_trajectory(80,1000+seed,False)
    val_np=generate_trajectory(20,2000+seed,False)
    ood_np=generate_trajectory(20,3000+seed,True)

    train=torch.tensor(train_np);val=torch.tensor(val_np);ood=torch.tensor(ood_np)
    P=primitive(train)
    mean=P.mean(dim=(0,1,2,3)).numpy();std=P.std(dim=(0,1,2,3)).numpy()
    state_std=torch.tensor([float(train[...,j].std()) for j in range(3)])

    torch.manual_seed(4000+seed)
    base=Solver((mean,std),"plain")
    base_state=base.face_net.state_dict()

    rows=[]
    models={}
    for mode in ["plain","soft","hard"]:
        model=Solver((mean,std),mode)
        model.face_net.load_state_dict(base_state)
        opt=torch.optim.Adam(model.parameters(),lr=5e-4)
        gen=torch.Generator().manual_seed(5000+seed)
        for _ in range(iterations):
            B=20
            inds=torch.randint(0,train.shape[0],(B,),generator=gen)
            ts=torch.randint(0,train.shape[1]-1,(B,),generator=gen)
            U=torch.stack([train[i,t] for i,t in zip(inds,ts)],0)
            target=torch.stack([train[i,t+1] for i,t in zip(inds,ts)],0)
            pred,pen=model.one_step(U)
            loss=(((pred-target)/state_std)**2).mean()+model.soft_weight*pen
            opt.zero_grad();loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(),1.0)
            opt.step()
        model.eval();models[mode]=model
        torch.save(model.state_dict(),out/f"seed{seed}_{mode}.pt")

    @torch.no_grad()
    def evaluate(model,data,split):
        U=data[:,0].clone();errs=[];minh=1e9;vc=0;vt=0;vmax=0.
        for t in range(1,data.shape[1]):
            Fx,Gy=model.fluxes(U)
            rx=entropy_residual(Fx,U,"x");ry=entropy_residual(Gy,U,"y")
            vc+=int((rx>1e-5).sum())+int((ry>1e-5).sum())
            vt+=rx.numel()+ry.numel()
            vmax=max(vmax,float(rx.max()),float(ry.max()))
            U=fv_step(U,Fx,Gy)
            minh=min(minh,float(U[...,0].min()))
            errs.append(float((((U-data[:,t])/state_std)**2).mean()))
        return dict(seed=seed,split=split,method=model.mode,
                    rollout_nrmse=float(np.sqrt(np.mean(errs))),min_h=minh,
                    entropy_violation_rate=vc/max(1,vt),max_entropy_residual=vmax)

    @torch.no_grad()
    def baseline(data,split):
        U=data[:,0].clone();errs=[];minh=1e9
        for t in range(1,data.shape[1]):
            U=fv_step(U,t_rusanov_x(U),t_rusanov_y(U))
            minh=min(minh,float(U[...,0].min()))
            errs.append(float((((U-data[:,t])/state_std)**2).mean()))
        return dict(seed=seed,split=split,method="rusanov",
                    rollout_nrmse=float(np.sqrt(np.mean(errs))),min_h=minh,
                    entropy_violation_rate=0.,max_entropy_residual=np.nan)

    for split,data in [("ID",val),("OOD",ood)]:
        for m in models.values(): rows.append(evaluate(m,data,split))
        rows.append(baseline(data,split))

    df=pd.DataFrame(rows);df.to_csv(out/f"seed{seed}_metrics.csv",index=False)
    np.savez(out/f"seed{seed}_stats.npz",mean=mean,std=std,state_std=state_std.numpy())
    return df


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--seed",type=int,required=True)
    p.add_argument("--outdir",default="/mnt/data/swe_2d_hcfl_runs")
    p.add_argument("--iterations",type=int,default=650)
    a=p.parse_args()
    print(train_seed(a.seed,a.outdir,a.iterations).to_string(index=False))
