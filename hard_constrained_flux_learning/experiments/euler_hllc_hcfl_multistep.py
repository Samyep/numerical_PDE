
import argparse, sys
from pathlib import Path
import numpy as np
import torch

sys.path.insert(0,"/mnt/data")
import euler_1d_hcfl as E
import euler_hllc_hcfl as H
import euler_1d_phase2_broad as B

torch.set_num_threads(2)

def finetune(seed,outdir="/mnt/data/euler_hllc_hcfl_multistep",iters=700,K=4):
    out=Path(outdir);out.mkdir(exist_ok=True)

    # Broad trajectory mixture.
    idata=torch.tensor(E.generate_trajectory(220,6000+seed,ood=False))
    odata=torch.tensor(E.generate_trajectory(260,7000+seed,ood=True))
    xdata=torch.tensor(B.generate_extreme_trajectory(100,8000+seed))
    train=torch.cat([idata,odata,xdata],0)

    st=np.load(f"/mnt/data/euler_hllc_hcfl/seed{seed}_broad_stats.npz")
    model=H.Solver((st["mean"],st["std"]))
    model.load_state_dict(torch.load(f"/mnt/data/euler_hllc_hcfl/seed{seed}_broad.pt",map_location="cpu"))
    model.train()

    state_std=torch.tensor(st["state_std"],dtype=torch.float32)
    opt=torch.optim.Adam(model.parameters(),lr=1.5e-4)
    gen=torch.Generator().manual_seed(15000+seed)

    for _ in range(iters):
        Bsz=36
        inds=torch.randint(0,train.shape[0],(Bsz,),generator=gen)
        ts=torch.randint(0,train.shape[1]-K,(Bsz,),generator=gen)
        U=torch.stack([train[i,t] for i,t in zip(inds,ts)],0)
        loss=0.0
        for k in range(1,K+1):
            target=torch.stack([train[i,t+k] for i,t in zip(inds,ts)],0)
            U=E.fv_step(U,model.flux(U))
            loss=loss+(((U-target)/state_std)**2).mean()/K

        opt.zero_grad();loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(),1.0)
        opt.step()

    model.eval()
    torch.save(model.state_dict(),out/f"seed{seed}_broad_multistep.pt")
    np.savez(out/f"seed{seed}_stats.npz",mean=st["mean"],std=st["std"],state_std=st["state_std"])
    return model

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--seed",type=int,required=True)
    p.add_argument("--iters",type=int,default=700)
    p.add_argument("--outdir",default="/mnt/data/euler_hllc_hcfl_multistep")
    a=p.parse_args()
    finetune(a.seed,a.outdir,a.iters)
    print("done",a.seed)
