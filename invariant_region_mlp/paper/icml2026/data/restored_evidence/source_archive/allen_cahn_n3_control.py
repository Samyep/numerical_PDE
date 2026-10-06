import numpy as np, math, json, sys, time
REF=.052802; T=.3; d=100

def g(x): return 1/(2+.4*np.dot(x,x))
def f(u): return u-u**3
class MLP:
 def __init__(self,M,mode,seed): self.M=M; self.mode=mode; self.rng=np.random.default_rng(seed); self.viols=[]
 def proj(self,u):
  self.viols.append(float((u<0) or (u>1)))
  if self.mode=='hard': return min(1.,max(0.,u))
  if self.mode=='clip': return float(np.clip(u,-1,1))
  return u
 def comp(self,t,x,n):
  if n==0:return 0.
  Mn=self.M**n; dt=T-t
  Z=self.rng.normal(size=(Mn,d)); XT=x[None,:]+math.sqrt(2*dt)*Z
  u=sum(g(xx) for xx in XT)/Mn
  for l in range(n):
   m=self.M**(n-l); b=0.
   for i in range(m):
    R=t+dt*self.rng.uniform(); Y=x+math.sqrt(2*(R-t))*self.rng.normal(size=d)
    ul=self.proj(self.comp(R,Y,l)); fl=f(ul)
    if l>0:
     um=self.proj(self.comp(R,Y,l-1)); fl-=f(um)
    b+=fl
   u+=dt*b/m
  return self.proj(u)
rows=[]
for M in [2,3,4]:
 for mode in ['raw','hard','clip']:
  vals=[]; vr=[]; sec=[]
  for r in range(100):
   s=time.time(); m=MLP(M,mode,10000+100*M+r); y=m.comp(0,np.zeros(d),3); vals.append(y); vr.append(np.mean(m.viols)); sec.append(time.time()-s)
  vals=np.array(vals); row={'n':3,'M':M,'mode':mode,'mean':vals.mean(),'std':vals.std(ddof=1),'mae':np.mean(abs(vals-REF)),'rel_mae':np.mean(abs(vals-REF))/REF,'viol':np.mean(vr),'seconds':np.mean(sec)}; print(row,flush=True); rows.append(row)
open('/mnt/data/allen_cahn_n3_control.json','w').write(json.dumps(rows,indent=2))
