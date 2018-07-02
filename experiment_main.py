from Models import *
import ModelFitter as mf
import PidTuner as pt
import numpy as np
from csv_rw import read_model_csv
import matplotlib.pyplot as plt
import os

# Generate test data, write to file
# lin_model_para=dict()
# lin_model_para['name']  = 'model_1'
# lin_model_para['A']     = np.matrix([[0,1],[-2,-3]])
# lin_model_para['B']     = np.matrix([[0],[1]])
# lin_model_para['C']     = np.matrix([[2,0]])
# lin_model_para['D']     = 0

# PID_PARA=dict()
# PID_PARA['P']=0.
# PID_PARA['I']=0.
# PID_PARA['D']=0.
# lin = LIN_Siso_Model(y_ref=y_ref_step,PID_PARA=PID_PARA,lin_model_para=lin_model_para) 
# lin.integrate(y0=[0,0,0],t0=0,tf=10,delta_t=0.01)

# Load test data

t,u,y,r = read_model_csv('D:\\Syno\\PID-Tuner\\ModelTesting\\model_1.csv')

# fit model

order = 2

print("Fit model with N = ", len(t) ," ...")
m_fitter = mf.ModelFitter(order) # input lag / output lag
m_fitter.fitModel(y, u,t)

# Forw. sim. of fitted modle using measured inputs:

#y_pred = np.zeros(len(y))
#y_pred[:y_lag] = y[:y_lag]
#for i in range(y_lag, len(y-y_lag)):
#	y_pred[i] = m_fitter.predict(u[i-u_lag:i+1], y_pred[i-y_lag:i])

# PID fitting:

#dt = np.mean(np.diff(t))

pid = pt.PidTuner(m_fitter)
pid.fitPid(y, r,u,[-100,100],t)

# Plotting:

#plt.plot(t, y)
#plt.plot(t, y_pred)
#plt.legend(['y', 'y_pred'])
#plt.title('u_lag = ' + str(u_lag) + ' y_lag = ' + str(y_lag)) 
#plt.show()


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