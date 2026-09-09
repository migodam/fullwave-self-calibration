"""Readable manuscript panels from locked exploratory records, no rerun/tuning."""
from pathlib import Path
import json
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
SRC=HERE.parents[1]/"delegated/a2_highdim"
sys.path.insert(0,str(SRC))
from common import make_model, make_scene

def main():
    rows=[json.loads(x) for x in (SRC/"records.jsonl").read_text().splitlines() if x.strip()]
    phi=make_model(32,"gaussian49").basis
    truth=phi@np.asarray(make_scene(801)["alpha_true"])
    methods=[("coherent_fixedpose","Coherent, fixed wrong pose"),("coherent_joint","Coherent, joint"),("intensity_joint","Matched intensity, joint"),("oracle","Known-pose reference")]
    panels=[("Truth",truth)]
    for method,label in methods:
        r=next(r for r in rows if r["phase"]=="test" and r["seed"]==801 and r["radius_index"]==0 and r["method"]==method)
        a=phi@np.asarray(r["alpha_est"])
        panels.append((label+f"\nRMSE {np.sqrt(np.mean((a-truth)**2)):.4f}",a))
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    fig,axs=plt.subplots(2,3,figsize=(9,6.5),layout="constrained")
    for ax,(label,a) in zip(axs.flat,panels):
        im=ax.imshow(a.reshape(32,32).T,origin="lower",extent=(-.5,.5,-.5,.5),vmin=truth.min(),vmax=truth.max(),cmap="viridis")
        ax.set(title=label,xlabel="x (m)",ylabel="y (m)")
        ax.set_xticks([-.5,0,.5]);ax.set_yticks([-.5,0,.5])
    axs[1,2].axis("off")
    axs[1,2].text(.05,.85,"Exploratory 2D test\nSeed 801; initial error 0.125 wavelengths\n\nDisplayed field: real contrast\n(reconstructed on the data grid).\n\nKnown pose is an oracle reference,\nnot a guaranteed performance bound.\n\nAll 12 test cases per method\ncontribute to the manuscript table.",va="top",fontsize=10)
    fig.colorbar(im,ax=list(axs.flat[:5]),shrink=.75,label=r"$\operatorname{Re}(\epsilon_r)-1$")
    fig.savefig(HERE/"figures/reconstruction_map_paper.png",dpi=200)
    plt.close(fig)

if __name__=="__main__": main()
