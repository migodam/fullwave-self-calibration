"""Plot actual campaign geometry without solving fields."""
from pathlib import Path
import sys,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
c=json.loads((ROOT/'runs/a2/imaging/frozen_config.json').read_text())['config']
v=VIE(Geometry(n=c['inverse_n'],n_tx=c['n_tx'],n_rx=c['n_rx'],aperture=c['aperture']),1.5e9)
fig,ax=plt.subplots(figsize=(8.2,5.7),layout='constrained');train=np.arange(len(v.rx))%2==0
ax.add_patch(Rectangle((-20,-20),40,40,facecolor='0.96',edgecolor='0.35',label='40 cm inversion region'))
ax.scatter(v.tx[:,0]*100,v.tx[:,1]*100,marker='^',s=70,label='12 transmitters (full circle)')
ax.scatter(v.rx[train,0]*100,v.rx[train,1]*100,s=22,label='24 fit receivers')
ax.scatter(v.rx[~train,0]*100,v.rx[~train,1]*100,marker='x',s=22,label='24 held receivers')
ax.text(0,0,'N128 integrated data\nN64 point inverse\n1.5 / 2.75 GHz',ha='center',va='center',fontsize=10)
ax.set(xlabel='x (cm)',ylabel='y (cm)',xlim=(-40,40),ylim=(-40,40),aspect='equal',title='Actual A2 synthetic acquisition')
ax.legend(loc='upper left',bbox_to_anchor=(1.01,1),fontsize=8);fig.savefig(ROOT/'figures/a2/acquisition.png',dpi=150)
