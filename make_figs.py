import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fractions import Fraction as Fr

INK, INK2, GRID = "#0b0b0b", "#52514e", "#dcdad4"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.family":"serif","font.serif":["DejaVu Serif"],"mathtext.fontset":"cm","font.size":10,"axes.edgecolor":INK2,
    "axes.labelcolor":INK,"xtick.color":INK2,"ytick.color":INK2,"axes.linewidth":0.6,
    "xtick.major.width":0.6,"ytick.major.width":0.6})

R=np.array([1.,2.,4.]); W=3.0; H=12/7; A=7/3
def S(t): return np.sum(1.0/(t+R))
def g(t): return W/S(t)-t
# interior branch, parametrized by t in (0, inf)
ts=np.concatenate([np.linspace(1e-6,1,400),np.linspace(1,400,4000)])
rho_i=np.array([g(t) for t in ts])
x_i=np.array([[(t+rh)/(W*(t+r)) for r in R] for t,rh in zip(ts,rho_i)])
pa_i=np.array([W*t/(t+rh) for t,rh in zip(ts,rho_i)])
pb_i=np.array([W/(t+rh) for t,rh in zip(ts,rho_i)])
lo=np.linspace(1.35,H,200); hi=np.linspace(A,2.55,200)

def style(ax):
    ax.grid(True,color=GRID,lw=0.5); ax.set_axisbelow(True)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    for v,lab in ((H,r"$H=12/7$"),(A,r"$A=7/3$")):
        ax.axvline(v,color=INK2,lw=0.6,ls=":")
        ax.text(v,1.01,lab,transform=ax.get_xaxis_transform(),ha="center",va="bottom",color=INK2,fontsize=9)
    ax.set_xlim(1.35,2.55); ax.set_xlabel(r"Capacity ratio $\rho=c_b/c_a$")

# ---- Figure 1: service shares
fig,ax=plt.subplots(figsize=(6.2,3.6))
cols=[S1,S2,S3]; lss=["-","--","-."]; names=[r"$r=1$",r"$r=2$",r"$r=4$"]
for j in range(3):
    xs=np.concatenate([lo,rho_i[::-1],hi])
    ys=np.concatenate([lo/(W*R[j]),x_i[::-1,j],np.full_like(hi,1/W)])
    order=np.argsort(xs); ax.plot(xs[order],ys[order],color=cols[j],ls=lss[j],lw=1.5,label=names[j])
# direct labels (text in ink, placed at left edge next to each line)
for j,yv in enumerate([1.35/(W*R[j]) for j in range(3)]):
    ax.text(1.36,yv+0.012,names[j],color=INK,fontsize=9,va="bottom")
style(ax); ax.set_ylabel("Share of first resource $x_i$")
ax.legend(frameon=False,loc="upper right",fontsize=9,handlelength=2.6,labelcolor=INK)
fig.tight_layout(); fig.savefig("/home/claude/unified/figs/fig_shares.pdf"); fig.savefig("/home/claude/unified/figs/fig_shares.png",dpi=160)

# ---- Figure 2: prices
fig,ax=plt.subplots(figsize=(6.2,3.4))
xs=np.concatenate([lo,rho_i[::-1],hi])
pa=np.concatenate([np.zeros_like(lo),pa_i[::-1],np.full_like(hi,W)])
pb=np.concatenate([W/lo,pb_i[::-1],np.zeros_like(hi)])
o=np.argsort(xs)
ax.plot(xs[o],pa[o],color=S1,ls="-",lw=1.5,label=r"$p_a$ (fixed resource)")
ax.plot(xs[o],pb[o],color=S2,ls="--",lw=1.5,label=r"$p_b$ (expanded resource)")
ax.text(2.40,2.80,r"$p_a$",color=INK,fontsize=9); ax.text(1.38,2.25,r"$p_b$",color=INK,fontsize=9)
style(ax); ax.set_ylabel("Resource price")
ax.legend(frameon=False,loc="center left",bbox_to_anchor=(0.0,0.40),fontsize=9,handlelength=2.6,labelcolor=INK)
fig.tight_layout(); fig.savefig("/home/claude/unified/figs/fig_prices.pdf"); fig.savefig("/home/claude/unified/figs/fig_prices.png",dpi=160)
print("figures written")
