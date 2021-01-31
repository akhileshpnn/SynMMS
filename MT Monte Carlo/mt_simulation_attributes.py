import numpy as np

def force_from_membrane(dh): 
    c1=0.06 # pN μm-1
    Fe=c1*dh
    return Fe 

def velocity_g(Fe,Tub): ## growth velocity
    c2=0.005 #μm s-1 μM-1
    vg=c2*Tub*np.exp(-Fe)
    return vg
   
def cat_freq_update(Tub):
    c3=0.04 #s-1
    c4=0.1 # μM-1
    kc = c3*np.exp(-c4*Tub)
    return kc

def cat_prob_update(l,kc, Rcell, dt):
    if l<=Rcell:
        factor=0
    elif l>Rcell:
        factor=1
    new_pc = 1-np.exp(-kc*dt) + factor*(1-np.exp(-kc*0.5*l*dt)) + factor*((1-np.exp(-kc*dt))*(1-np.exp(-kc*0.5*l*dt)))
    return new_pc

def attractive_force(x,h,dh,deforming_len): ## attractive force that depends on distance and deformation
    if deforming_len<=x:
        Ia=0
    elif deforming_len>x:
        Ia=dh*(1-x/deforming_len)/h
    return abs(Ia)

def finding_attracting_mt(mt_id,Force_of_attraction, Nmtmax): # returns the mt that attracts the given mt_id
    force=[]
    for q in range(len(Force_of_attraction)):        
        if mt_id==Force_of_attraction[q][1]:
            force.append(Force_of_attraction[q][2])
        else:
            force.append(0)
    if force.count(0)==len(force):       
        return Nmtmax*10
    else:
        attracting_mt=Force_of_attraction[np.argmax(force)][0]  ## works  -- double checked
        return attracting_mt

def finding_attracting_mt_with_largest_length(mt_id,Force_of_attraction,Lmt): # returns the mt with largest length that attracts the given mt_id
    force = []
    attracting_mts = []
    for q in range(len(Force_of_attraction)):        
        if mt_id==Force_of_attraction[q][1]:
            attracting_mts.append(Force_of_attraction[q][0])
            force.append(Force_of_attraction[q][2])
    Len_attracting_mts = []
    for i in range(len(attracting_mts)):
        Len_attracting_mts.append(Lmt[attracting_mts[i]])
    mt_id_with_largest_defor = np.argmax(Len_attracting_mts)
    attracting_mt_with_largest_length = attracting_mts[mt_id_with_largest_defor]
    force_by_domi = force[mt_id_with_largest_defor]
    return attracting_mt_with_largest_length,force_by_domi
        
## to find the for ce of attraction by one mt on another

def finding_force_of_attraction(attracting_mt_id,attracted_mt_id,Force_of_attraction,Deforming_MT_id):
    for q in range(len(Force_of_attraction)):
        if Deforming_MT_id[attracting_mt_id] == Force_of_attraction[q][0] and Deforming_MT_id[attracted_mt_id] == Force_of_attraction[q][1]:
            ff = Force_of_attraction[q][2]
    return ff
    
def Attracted_MT(deforming_mt_id,min_radius,prev_position, Nmtmax):
    
    neighbours=[]
    
    x0=prev_position[deforming_mt_id]
    for i in range(Nmtmax):
        x=prev_position[i]
        dist=abs(x0-x)
        if x0<10+np.pi:
            dist_around=abs(x0+2*np.pi-x)
        else:
            dist_around=abs((x0-2*np.pi)-x)
        min_dist=min(dist,dist_around)
        if x0==x:
            pass
        elif min_dist<min_radius:
            neighbours.append(i) 

    return neighbours
            


def checking_for_element(element,array):
    element = int(element)
    out = False
    for i in range(len(array)):
        if int(array[i])==element:
            out = True
            break
        else:
            continue
    return out
## finding the index of an element in an array
def Finding_index(element,array):
    
    for i in range(len(array)):
        if array[i]==element:
            return i
            break

## pairwise sorting of deforming mts

def Pairwise_sorting(Deforming_mts):
    Pairwise = []    
    l = len(Deforming_mts) 
    if l<=1:
        out = False
    elif l>1:
        m = 0
        while m<l:
            for j in range(m+1,l):
                Pairwise.append([Deforming_mts[m],Deforming_mts[j]])
            m = m+1
        out = Pairwise
    return out


fraction = 0.9
def MT_position_update(deforming_mt,Prev_position,Attracted_neighbour):

    New_position = [Prev_position[i] for i in range(len(Prev_position))]
    
    x0=Prev_position[deforming_mt]
    for i in range(len(Attracted_neighbour)):
        attracted_mt_indx = Attracted_neighbour[i]
        x=Prev_position[attracted_mt_indx]
        dist=abs(x0-x)
        if x0<10+np.pi:
            dist_around=abs(x0+2*np.pi-x)
        else:
            dist_around=abs((x0-2*np.pi)-x)
        min_dist=min(dist,dist_around)
        
        new_distance = min_dist*fraction
        if x0==x:
            pass
        if x>x0 and min_dist==dist_around:
            New_position[attracted_mt_indx] = (2*np.pi+New_position[deforming_mt])-new_distance
            if New_position[attracted_mt_indx]>10+2*np.pi:
                 New_position[attracted_mt_indx]=abs(New_position[attracted_mt_indx]-2*np.pi)
        elif x>x0 and min_dist==dist:
            New_position[attracted_mt_indx] = New_position[deforming_mt]+new_distance 
#            if New_position[attracted_mt_indx]<10:
#                 New_position[attracted_mt_indx]=abs(New_position[attracted_mt_indx]+2*np.pi)
        elif x<x0 and min_dist==dist_around:
            New_position[attracted_mt_indx] = New_position[deforming_mt]-2*np.pi+new_distance 
            if New_position[attracted_mt_indx]<10:
                 New_position[attracted_mt_indx]=abs(New_position[attracted_mt_indx]+2*np.pi)
        elif x<x0 and min_dist==dist:
            New_position[attracted_mt_indx] = New_position[deforming_mt]-new_distance 
#            if New_position[attracted_mt_indx]<10:
#                 New_position[attracted_mt_indx]=abs(New_position[attracted_mt_indx]+2*np.pi)
    return New_position

    

def updating_neighbours(Attracted_neighbours, LLmt, Rcell, time_point):        
#        # Removing neighbours that didnt touch the membrane

        for c in range(len(Attracted_neighbours)):
            Index_for_deleting_non_deforming = []
            if len(Attracted_neighbours[c])==0:
                continue
            else:
                for d in range(len(Attracted_neighbours[c])):
                    if LLmt[:,time_point][Attracted_neighbours[c][d]]<=Rcell:
                        Index_for_deleting_non_deforming.append(Finding_index(Attracted_neighbours[c][d],Attracted_neighbours[c]))                        
                    else:
                        continue
            Attracted_neighbours[c] = np.delete(Attracted_neighbours[c],Index_for_deleting_non_deforming)  
        return Attracted_neighbours

def sorting_mts(Attracted_neighbours, Deforming_MT_id, Deforming_len):    
    Deforming_len_sortedIndxs = np.argsort(Deforming_len)
    Deforming_MT_id=np.array(Deforming_MT_id)
    Attracted_neighbours=np.array(Attracted_neighbours)
    if len(Deforming_len_sortedIndxs)!=0:
        Deforming_MT_id=Deforming_MT_id[Deforming_len_sortedIndxs]
        Attracted_neighbours=Attracted_neighbours[Deforming_len_sortedIndxs]

    
    return Deforming_MT_id, Attracted_neighbours
