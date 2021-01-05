import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set up plotting
rc={'image.cmap': 'YlGnBu', 'font.size': 20, 'axes.labelsize': 20, 'legend.fontsize': 10.0, 
    'axes.titlesize': 30, 'xtick.labelsize': 20, 'ytick.labelsize': 20, 'lines.linewidth': 5.0,
    'fontname':'Arial'}
sns.set_context(rc=rc)


def plot_bundle(t,MT_positions,MT_heights,xrange,Rcell,label=None): #xrange is the artificial spread to MTs to make the bundles visible in snapshot
    
    delta_distance=1e-2
    
    mt_positions=MT_positions[t]
    mt_heights=MT_heights[:,t]
    neighbour_sep=mt_positions[1:]-mt_positions[:-1]
    num_of_bundles=len(np.argwhere(neighbour_sep>delta_distance))+1
    bundle_positions=mt_positions[np.argwhere(neighbour_sep>delta_distance)]
    bundle_size=np.argwhere(neighbour_sep>delta_distance)+1
    if num_of_bundles>1:
        bundle_positions=np.append(bundle_positions,mt_positions[-1])
        bundle_size=np.append(bundle_size,len(mt_positions))
        bundle_size=np.pad(bundle_size,pad_width=(1, 0), mode='constant')
    plt.figure()
    if t==1:
        plt.bar(mt_positions, mt_heights, color='goldenrod', width = 0.02)    
    
    else:
        for i in range(num_of_bundles):    
            bundle_pos = np.linspace(bundle_positions[i]-xrange, bundle_positions[i]+xrange,bundle_size[i+1]-bundle_size[i])
            bundle_heights = mt_heights[bundle_size[i]:bundle_size[i+1]]
            plt.bar(bundle_pos, bundle_heights, color='goldenrod', width = 0.02)

    plt.axhline(y=Rcell,color='b',ls='--',lw=3.0)
    plt.xticks([10,10+np.pi,10+2*np.pi],[0, '$\pi$', '2$\pi$'])
    plt.xlabel('MT-position')
    if label==True:
        
        plt.ylabel('MT-length($\mu m$)')     
    else:
        plt.yticks([])
    plt.ylim(0,50)
    plt.xlim(9,17)
#    plt.show()
