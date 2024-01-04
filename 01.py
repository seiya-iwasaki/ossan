import numpy as np
import matplotlib.pyplot as plt
import csv
from csv import writer
import math

base_depth = 1.5
fl = 0.6
floor_hight = 3
floor_weight_force = 0 #ini 6246
sq_unit_floor_weight_force = 13
unit_base_weight_force = 0 #ini
unit_floor_weight_force = 0 #ini
reduction_rate =  0.4
opening_bottom = 0.0
opening_top = 2.2
air_pool_hight = floor_hight-opening_top
effective_float_weight_ratio = 0.42
rho_g = 9.8
finess = 0.01
open_cell = [12.1,16.125]
cell_length = [6,8]

index_a = [3.0,2.0,1.5]
index_mu = [0.4,0.5,0.6] 

def truncate(number, digits) -> float:
    # Improve accuracy with floating point operations, to avoid truncate(16.4, 2) = 16.39 or truncate(-1.13, 2) = -1.12
    nbDecimals = len(str(number).split('.')[1]) 
    if nbDecimals <= digits:
        return number
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper


def baloon(x_lst,S,W,H):
    y_lst = []
    dx = x_lst[1]-x_lst[0]
    for i, x in enumerate(x_lst):
        if x > H+floor_hight/2:
            y_lst.append(y_lst[-1])
        elif x > H:
            top_float = effective_float_weight_ratio*unit_floor_weight_force*dx
            y_lst.append(y_lst[-1]-top_float)
        elif x<=0 : 
            y_lst.append(W)
        else:
            y = W-(rho_g*S*(x+base_depth))
            y_lst.append(y)

    for i, y in  enumerate(y_lst):
        if y < 0 : y_lst[i] = 0
        
    return np.array(y_lst)

def partial_sunk(x_lst,S,W,H):
    y_lst = []
    dx = x_lst[1]-x_lst[0]
    for i, x in enumerate(x_lst):
        base_float = rho_g*S*(fl+base_depth)
        if x > H+floor_hight/2:
            y_lst.append(y_lst[-1])
        elif x > H:
            top_float = effective_float_weight_ratio*unit_floor_weight_force*dx
            y_lst.append(y_lst[-1]-top_float)
        elif x<=0 : 
            y_lst.append(W)
        elif x < fl :
            base_float = rho_g*S*(x+base_depth)
            y = W-base_float
            y_lst.append(y)
        elif x < fl+(floor_hight/2):
            base_float += effective_float_weight_ratio*unit_base_weight_force*(x-fl)
            y = W-base_float
            y_lst.append(y)
        elif x < fl+opening_top:
            base_float += effective_float_weight_ratio*unit_base_weight_force*(floor_hight/2)
            first_float = effective_float_weight_ratio*unit_floor_weight_force*(x-floor_hight/2-fl)
            y = W-base_float-first_float
            y_lst.append(y)
        else:
            f_number = (x-(fl+opening_top)) // floor_hight
            f_x = (x-(fl+opening_top)) % floor_hight

            base_float += effective_float_weight_ratio*unit_base_weight_force*(floor_hight/2)
            first_float = effective_float_weight_ratio*unit_floor_weight_force*(opening_top-floor_hight/2)

            upper_floor_float = f_number*((effective_float_weight_ratio*unit_floor_weight_force*opening_top)+(rho_g*S*air_pool_hight))
            y = W - base_float - first_float - upper_floor_float

            if f_x<=air_pool_hight+opening_bottom:
                y -= rho_g*S*f_x
            else:
                y -= (effective_float_weight_ratio*unit_floor_weight_force*(f_x-air_pool_hight))+(rho_g*S*air_pool_hight)
            y_lst.append(y)

    for i, y in  enumerate(y_lst):
        if y < 0 : y_lst[i] = 0
    
    return np.array(y_lst)

def all_sunk(x_lst,S,W,H):
    y_lst = []
    dx = x_lst[1]-x_lst[0]
    for i, x in enumerate(x_lst):
        base_float = rho_g*S*(fl+base_depth)
        if x > H+floor_hight/2:
            y_lst.append(y_lst[-1])
        elif x > H:
            top_float = effective_float_weight_ratio*unit_floor_weight_force*dx
            y_lst.append(y_lst[-1]-top_float)
        elif x<=0 : 
            y_lst.append(W)
        elif x < fl :
            base_float = rho_g*S*(x+base_depth)
            y = W-base_float
            y_lst.append(y)
        elif x < fl+(floor_hight/2):
            base_float += effective_float_weight_ratio*unit_base_weight_force*(x-fl)
            y = W-base_float
            y_lst.append(y)
        elif x < fl+opening_top:
            base_float += effective_float_weight_ratio*unit_base_weight_force*(floor_hight/2)
            first_float = effective_float_weight_ratio*unit_floor_weight_force*(x-floor_hight/2-fl)
            y = W-base_float-first_float
            y_lst.append(y)
        else:
            f_number = (x-(fl+opening_top)) // floor_hight
            f_x = (x-(fl+opening_top)) % floor_hight
            
            base_float += effective_float_weight_ratio*unit_base_weight_force*(floor_hight/2)
            first_float = effective_float_weight_ratio*unit_floor_weight_force*(opening_top-floor_hight/2)

            upper_floor_float = f_number*(effective_float_weight_ratio*unit_floor_weight_force*floor_hight)

            y = W - base_float - first_float - upper_floor_float

            if f_x<=air_pool_hight+opening_bottom:
                y -= rho_g*S*f_x
            else:
                y -= (effective_float_weight_ratio*unit_floor_weight_force*f_x)
            y_lst.append(y)

    for i, y in  enumerate(y_lst):
        if y < 0 : y_lst[i] = 0

    return np.array(y_lst)


def tsunami_force_func(x_lst, a, B, H):
    y_lst = []
    for i, x in enumerate(x_lst):
        if x*a < H:
            y_lst.append((1/2)*rho_g*B*a*a*(x*x))
        else:
            y_lst.append((1/2)*rho_g*B*(2*a*x*H-(H*H)))

    return np.array(y_lst)


with open('conditions.csv', 'r') as cfd, open('points.csv', 'w+') as pfd:

    reader = csv.reader(cfd)
    writer = writer(pfd)

    counter  = 0
    for conditions in reader:

        figure, axis = plt.subplots(3, 2)
        
        print(conditions)
        for i in range(4):
            conditions[i] = int(conditions[i])
        
        base_depth = float(conditions[4])

        cell_num = [conditions[2] * (conditions[0] / cell_length[0]),conditions[2] * (conditions[1] / cell_length[1])]

        S = conditions[0]*conditions[1] #床面積

        floor_num = conditions[2]

        floor_weight_force = S*sq_unit_floor_weight_force

        weight_reduction_rate = 1 - 2*(cell_num[0]*open_cell[0]+cell_num[1]*open_cell[1])/(S*(floor_num+1)+2*(cell_num[0]*open_cell[0]+cell_num[1]*open_cell[1]))

        unit_floor_weight_force = floor_weight_force/floor_hight
        unit_base_weight_force = floor_weight_force/(base_depth+fl+floor_hight/2)
        W = floor_weight_force*(conditions[2]+1) 

        if(conditions[3]==1): #skelton case
            unit_floor_weight_force *= weight_reduction_rate
            unit_base_weight_force *= weight_reduction_rate
            W *= weight_reduction_rate

        H = floor_hight*conditions[2]+fl

        points = []

        for i in range(2): # X or Y

            B = conditions[i] #0,1->X,Y graph width

            x = np.arange(0, 10, finess)

            effective_weight_baloon = baloon(x,S,W,H)
            effective_weight_normal = partial_sunk(x,S,W,H)
            effective_weight_sunk = all_sunk(x,S,W,H)

            for k in range(3):
                for a in index_a: #culcurating tsunami force

                    tsunami_force = tsunami_force_func(x, a, B+0.5, H) 

                    if(conditions[3]==1): #skelton case

                        all_area = (B+0.5)*H
                        open_area = open_cell[i] *  cell_num[i]

                        reduction_rate = (all_area - open_area)/all_area

                        tsunami_force = tsunami_force*reduction_rate
                

                     #print(len(tsunami_force))
                    axis[k, i].plot(x, tsunami_force,label=rf'a = {a}')            
                
                    for mu in index_mu:
                    
                        if k  == 0:
                            idx = np.argwhere(np.diff(np.sign(tsunami_force - mu*effective_weight_baloon))).flatten()
                        elif k == 1:
                            idx = np.argwhere(np.diff(np.sign(tsunami_force - mu*effective_weight_normal))).flatten()
                        elif k==2:
                            idx = np.argwhere(np.diff(np.sign(tsunami_force - mu*effective_weight_sunk))).flatten()
                        
                        if len(idx) != 0:
                            #print(idx)
                            points.append(truncate(x[idx[0]],2))
                            #points.append(truncate(tsunami_force[idx[0]],0))
                            #print(truncate(x[idx[0]],2))
                        else:
                            points.append("None")
                            #points.append("None")
                            #print("None")
                        
                
            for A, mu in enumerate(index_mu):
                axis[0, i].plot(x, mu*effective_weight_baloon,label=rf'$\mu = {mu} $')
                axis[1, i].plot(x, mu*effective_weight_normal,label=rf'$\mu = {mu} $')
                axis[2, i].plot(x, mu*effective_weight_sunk,label=rf'$\mu = {mu} $')

        axis[0, 1].legend(loc='upper left',bbox_to_anchor=(1,1))
        axis[0,0].set_title("X")
        axis[0,1].set_title("Y")
        axis[2,0].set_xlabel("tsunami hight [m]")
        axis[2,1].set_xlabel("tsunami hight [m]")
        axis[0,0].set_ylabel("condition (a)\n force [kN]")
        axis[1,0].set_ylabel("condition (b)\n force [kN]")
        axis[2,0].set_ylabel("condition (c)\n force [kN]")
        for i in range(3):
            for j in range(2):
                axis[i,j].set_ylim(0,60000)
        
        writer.writerow(points)

    
    pfd.close()

    plt.show()

pfd.close()
