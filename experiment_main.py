from Models import *
import ModelFitter as mf
import numpy as np

lin_model_para=dict()
lin_model_para['name']  = 'model_1'
lin_model_para['A']     = np.matrix([[0,1],[-2,-3]])
lin_model_para['B']     = np.matrix([[0],[1]])
lin_model_para['C']     = np.matrix([[2,0]])
lin_model_para['D']     = 0

PID_PARA=dict()
PID_PARA['P']=10
PID_PARA['I']=5
PID_PARA['D']=1
lin = LIN_Siso_Model(y_ref=y_ref_step,PID_PARA=PID_PARA,lin_model_para=lin_model_para) 
lin.integrate(y0=[0,0,0],t0=0,tf=10,delta_t=0.001)

# lin_model_para=dict()
# lin_model_para['name']  = 'model_2'
# lin_model_para['A']     = np.matrix([[0,1,0],[0,0,1],[-1,-1,-1]])
# lin_model_para['B']     = np.matrix([[0],[0],[1]])
# lin_model_para['C']     = np.matrix([[1,1,1]])
# lin_model_para['D']     = 0
    

# PID_PARA=dict()
# PID_PARA['P']=10
# PID_PARA['I']=5
# PID_PARA['D']=1
# lin = LIN_Siso_Model(y_ref=y_ref_step,PID_PARA=PID_PARA,lin_model_para=lin_model_para) 
# lin.integrate(y0=[1,1,1,0],t0=0,tf=5,delta_t=0.001)
