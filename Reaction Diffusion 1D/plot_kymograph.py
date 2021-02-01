import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

# Set up plotting
rc={'image.cmap': 'YlGnBu', 'font.size': 20, 'axes.labelsize': 20, 'legend.fontsize': 10.0, 
    'axes.titlesize': 20, 'xtick.labelsize': 20, 'ytick.labelsize': 20, 'lines.linewidth': 3.0}
sns.set_context(rc=rc)

input_folder = os.path.abspath(os.getcwd())+'\\saved data5\\'

t_total = int(1.0e8) 

U1_array = [];U2_array = []
Light_array = []

for i in range(0,t_total,200000):
    print(i, end='\r')
    u1 = np.load(os.path.join(input_folder,'U1_loop_'+str(i)+'.npy'))
    u2 = np.load(os.path.join(input_folder,'U2_loop_'+str(i)+'.npy'))
    light = np.load(os.path.join(input_folder,'Light cue_loop_'+str(i)+'.npy'))
    
    U1_array.append(u1)
    U2_array.append(u2)
    Light_array.append(light)

U1_array=np.array(U1_array);
U2_array=np.array(U2_array);
Light_array=np.array(Light_array);

plt.figure()
plt.imshow(np.transpose(U1_array),cmap='viridis',aspect='auto',origin='lower')
if np.max(Light_array)>1.0:
    Masked_pheromone=Light_array.copy()
    Masked_pheromone[Light_array<=1.2]=np.nan
    plt.imshow(np.transpose(Masked_pheromone),cmap='gray',aspect='auto',origin='lower',alpha=0.4)

plt.yticks(ticks=[0,500,1000],labels=['0','180','360'])
plt.ylabel('Angle($^\circ$)')
plt.xlabel('Time(arb.units)')
if t_total==int(1.0e8):
    plt.xticks(ticks=[0,100,200,300,400,500],labels=['0','2000','4000','6000','8000','10000']) 
elif t_total==int(2.0e8):
    plt.xticks(ticks=[0,100,200,300,400,500,600,700,800,900,1000],labels=['0','2000','4000','6000','8000','10000','12000','14000','16000','18000','20000'])
plt.show()
