
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch import nn

sys.path.insert(0, "/mnt/data")
import euler_1d_hcfl as E

torch.set_num_threads(2)

def hllc_pair(UL, UR):
    rhoL, rhoR = UL[...,0].clamp_min(1e-10), UR[...,0].clamp_min(1e-10)
    mL, mR = UL[...,1], UR[...,1]
    EL, ER = UL[...,2], UR[...,2]
    uL, uR = mL/rhoL, mR/rhoR
    pL, pR = E.t_pressure(UL).clamp_min(1e-10), E.t_pressure(UR).clamp_min(1e-10)
    cL, cR = torch.sqrt(E.GAMMA*pL/rhoL), torch.sqrt(E.GAMMA*pR/rhoR)

    SL = torch.minimum(uL-cL, uR-cR)
    SR = torch.maximum(uL+cL, uR+cR)
    den = rhoL*(SL-uL) - rhoR*(SR-uR)
    SM = (pR-pL + rhoL*uL*(SL-uL) - rhoR*uR*(SR-uR)) / (den + 1e-14)

    rhoSL = rhoL*(SL-uL)/(SL-SM + 1e-14)
    rhoSR = rhoR*(SR-uR)/(SR-SM + 1e-14)
    ESL = rhoSL*(EL/rhoL + (SM-uL)*(SM + pL/(rhoL*(SL-uL) + 1e-14)))
    ESR = rhoSR*(ER/rhoR + (SM-uR)*(SM + pR/(rhoR*(SR-uR) + 1e-14)))
    USL = torch.stack([rhoSL, rhoSL*SM, ESL], dim=-1)
    USR = torch.stack([rhoSR, rhoSR*SM, ESR], dim=-1)

    fL, fR = E.t_flux(UL), E.t_flux(UR)
    FSL = fL + SL[...,None]*(USL-UL)
    FSR = fR + SR[...,None]*(USR-UR)

    return torch.where((SL>=0)[...,None], fL,
           torch.where((SM>=0)[...,None], FSL,
           torch.where((SR>0)[...,None], FSR, fR)))

def t_hllc(U):
    return hllc_pair(U, torch.roll(U,-1,dims=-2))

class HLLCLearnedFlux(nn.Module):
    def __init__(self, mean, std, width=72):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(15,width), nn.Tanh(),
            nn.Linear(width,width), nn.Tanh(),
            nn.Linear(width,3),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean,dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std,dtype=torch.float32))
        self.register_buffer("scale", torch.tensor([0.6,1.2,2.5],dtype=torch.float32))

    def forward(self,U):
        P=E.primitive(U)
        feats=torch.cat([(torch.roll(P,s,dims=-2)-self.mean)/self.std
                         for s in [2,1,0,-1,-2]], dim=-1)
        corr=torch.tanh(self.net(feats))
        PL,PR=P,torch.roll(P,-1,dims=-2)
        jump=torch.sqrt((((PR-PL)/self.std)**2).sum(dim=-1)+1e-12)
        return t_hllc(U) + 0.18*jump[...,None]*corr*self.scale

class Solver(nn.Module):
    def __init__(self, stats):
        super().__init__()
        self.flux_net=HLLCLearnedFlux(*stats)

    def flux(self,U):
        return E.hard_entropy_projection(self.flux_net(U),U)

    def one_step(self,U):
        return E.fv_step(U,self.flux(U))

def train_seed(seed,outdir="/mnt/data/euler_hllc_hcfl",iters=1100,broad=False):
    out=Path(outdir); out.mkdir(exist_ok=True)

    if broad:
        import euler_1d_phase2_broad as B
        idata=torch.tensor(E.generate_trajectory(220,6000+seed,ood=False))
        odata=torch.tensor(E.generate_trajectory(260,7000+seed,ood=True))
        xdata=torch.tensor(B.generate_extreme_trajectory(100,8000+seed))
        train=torch.cat([idata,odata,xdata],dim=0)
    else:
        train=torch.tensor(E.generate_trajectory(420,1000+seed,ood=False))

    P=E.primitive(train)
    mean=P.mean(dim=(0,1,2)).numpy()
    std=P.std(dim=(0,1,2)).numpy()
    state_std=torch.tensor([float(train[...,j].std()) for j in range(3)])

    torch.manual_seed(12000+seed)
    model=Solver((mean,std))
    opt=torch.optim.Adam(model.parameters(),lr=3e-4)
    gen=torch.Generator().manual_seed(13000+seed)

    for _ in range(iters):
        B=56
        inds=torch.randint(0,train.shape[0],(B,),generator=gen)
        ts=torch.randint(0,train.shape[1]-1,(B,),generator=gen)
        U=torch.stack([train[i,t] for i,t in zip(inds,ts)],0)
        target=torch.stack([train[i,t+1] for i,t in zip(inds,ts)],0)
        pred=model.one_step(U)
        loss=(((pred-target)/state_std)**2).mean()
        opt.zero_grad(); loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(),1.0)
        opt.step()

    model.eval()
    tag="broad" if broad else "id"
    torch.save(model.state_dict(),out/f"seed{seed}_{tag}.pt")
    np.savez(out/f"seed{seed}_{tag}_stats.npz",mean=mean,std=std,state_std=state_std.numpy())
    return model,(mean,std),state_std

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--seed",type=int,required=True)
    p.add_argument("--outdir",default="/mnt/data/euler_hllc_hcfl")
    p.add_argument("--iters",type=int,default=1100)
    p.add_argument("--broad",action="store_true")
    a=p.parse_args()
    train_seed(a.seed,a.outdir,a.iters,a.broad)
    print("done",a.seed,a.broad)
