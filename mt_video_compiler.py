# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 14:41:40 2020

@author: nandan
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import animation


folder_read = os.path.abspath(os.getcwd())+'\\saved data\\'
Nmtmax=20
Rcell = 25


dt = 1
N = 2000

X_position = [];LLmt = []
for i in range(0,N,dt):
    print(i,end='\r')
    position = np.load(os.path.join(folder_read+'\\','Mts('+str(Nmtmax)+')_position_loop_'+str(i)+'.npy'))
    height = np.load(os.path.join(folder_read+'\\','Mts('+str(Nmtmax)+')_height_loop_'+str(i)+'.npy'))
    X_position.append(position)
    LLmt.append(height)


fig=plt.figure()
plt.ylim(0,50)
plt.xlim(min(X_position[0])-1,max(X_position[0])+1)
plt.axhline(y=Rcell,color='b',ls='--',lw=3.0)
plt.xticks(ticks=[10,10+2*np.pi],labels=[0,'2$\pi$'])
plt.xlabel('MT position')
plt.ylabel('MT length($\mu m$)')

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

X=X_position[0]
barcollection = plt.bar(X,LLmt[0],bottom=0,width=0.01,color='k')

def animate(i):
    y=LLmt[i+1]
    for j, b in enumerate(barcollection):
        b.set_height(y[j])
        b.set_x(X_position[i+1][j])
        b.set_y(1)
        b.set_color('goldenrod')

    plt.title(str(i))
    print(i, end="\r")

anim=animation.FuncAnimation(fig,animate,repeat=True,blit=False,frames=len(X_position)-1)
anim.save('MTs('+str(Nmtmax)+').mov', writer=writer)
plt.show()