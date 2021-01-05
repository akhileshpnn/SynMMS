import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

# Set up plotting
rc={'image.cmap': 'YlGnBu', 'font.size': 20, 'axes.labelsize': 20, 'legend.fontsize': 10.0, 
    'axes.titlesize': 20, 'xtick.labelsize': 20, 'ytick.labelsize': 20, 'lines.linewidth': 3.0}
sns.set_context(rc=rc)



input_folder = os.path.abspath(os.getcwd())+'\\saved data4\\'

t_total = int(1.0e8) 

U1_array = [];U2_array = []
Pheromone_array = []

for i in range(0,t_total,200000):
    print(i, end='\r')
    u1 = np.load(os.path.join(input_folder,'U1_loop_'+str(i)+'.npy'))
    u2 = np.load(os.path.join(input_folder,'U2_loop_'+str(i)+'.npy'))
    pheromone = np.load(os.path.join(input_folder,'Light cue_loop_'+str(i)+'.npy'))
    
    U1_array.append(u1)
    U2_array.append(u2)
    Pheromone_array.append(pheromone)

U1_array=np.array(U1_array);
U2_array=np.array(U2_array);
Pheromone_array=np.array(Pheromone_array);

plt.figure()
plt.imshow(np.transpose(U1_array),cmap='viridis',aspect='auto',origin='lower')
if np.max(Pheromone_array)>1.0:
    Masked_pheromone=Pheromone_array.copy()
    Masked_pheromone[Pheromone_array<=1.2]=np.nan
    plt.imshow(np.transpose(Masked_pheromone),cmap='gray',aspect='auto',origin='lower',alpha=0.4)

plt.yticks(ticks=[0,500,1000],labels=['0','180','360'])
plt.ylabel('Angle($^\circ$)')
plt.xticks(ticks=[0,100,200,300,400,500],labels=['0','2000','4000','6000','8000','10000'])
plt.xlabel('Time(arb.units)')
plt.show()
