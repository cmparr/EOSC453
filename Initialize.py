import numpy as np 

global Flux_in_4, Flux_out_4, Flux_in_9, Flux_out_9, M0_4, M0_9

## 4 BOX MODEL
# initial masses
M0_4 = np.array([725, 725, 110, 60]) # 4 box model M0 values 
# fluxes between boxes
F12_4 = 90
F21_4 = 90
F15_4 = 110
F51_4 = 55
F57_4 = 55
F71_4 = 55 #detritus decomposition
F72_4 =  0

Flux_in_4 = np.array([[0, F21_4, F51_4, F71_4], \
                    [F12_4, 0, 0, F72_4], \
                    [F15_4, 0, 0, 0],\
                    [0, 0, F57_4, 0]])
Flux_out_4 = Flux_in_4.T
Mass_Flux_total_4 = sum(Flux_in_4 - Flux_out_4)
if sum(Mass_Flux_total_4) == 0:
    print('4 Box Net flux is 0')


## 9 BOX MODEL
# initial masses
M0_9 = np.array([725, 725, 3, 37675, 110, 450, 60, 1350, 160])
# M0_9 = np.array([725, 725, 3, 37+625, 110, 450, 60, 1350, 160]) # TODO why is M4 37+625?
# fluxes between boxes
F21 = 90
F51= 55
F61 = 0 #1000 # deforestation
F71 = 50 # detritus decomposition
F81 = 3
F91 = 1
F12 = 89
F32 = 36
F42 = 42
F72 = 1
F23 = 40
F24 = 38
F34 = 4 # detritus
F15 = 110
F56 = 15
F57 = 40
F67 = 15
F78 = 3
F79 = 1

# Mostly zeros, start there
Flux_out_9 = np.zeros((9,9))
# These locations should all correspond with the coefficients on the flux value, minus one for python
Flux_out_9[1, 0] = F21; Flux_out_9[4, 0] = F51; Flux_out_9[5, 0] = F61; Flux_out_9[6, 0] = F71; Flux_out_9[7, 0] = F81; Flux_out_9[8, 0]= F91
Flux_out_9[0, 1] = F12; Flux_out_9[2, 1] = F32; Flux_out_9[3, 1] = F42; Flux_out_9[6, 1] = F72; Flux_out_9[1, 2] = F23; Flux_out_9[2, 3] = F34
Flux_out_9[0, 4] = F15; Flux_out_9[4, 5] = F56; Flux_out_9[4, 6] = F57; Flux_out_9[5, 6] = F67; Flux_out_9[6, 7] = F78; Flux_out_9[6, 8] = F79
Flux_out_9[1, 3] = F24

Flux_in_9 = Flux_out_9.T
# check for zero
Mass_Flux_total_9 = sum(Flux_in_9 - Flux_out_9)
if sum(Mass_Flux_total_9) == 0:
    print('9 Box Net flux is 0')


# The timespan to integrate over
t_start = 1500; t_end= 2200; n = 200 # TODO want 1800-2200 time interval
t = np.linspace(t_start, t_end, n) # some time span