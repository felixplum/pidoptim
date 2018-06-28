#from casadi import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import csv

from csv_rw import read_model_csv
import os

#cwd = os.path.dirname(__file__)
#print(cwd)
#t,u,y = read_model_csv('D:\\Syno\\PID-Tuner\\ModelTesting\\model_1.csv')



class Siso_Model(object):
    def __init__(self,y_ref,PID_PARA,name):
        self.name=name
        self.KP=PID_PARA['P']
        self.KI=PID_PARA['I']
        self.KD=PID_PARA['D']
        self.u = 0
        self.y_ref=y_ref
        self.e_last = 0
        self.u_min = -5000
        self.u_max = 5000
        
      
        
        
    def integrate(self,y0,t0,tf,delta_t):
        self.dt=delta_t
        ode_inte = ode(self.f).set_integrator('vode', method='adams', with_jacobian=False,max_step=delta_t)
        ode_inte.set_initial_value(y0, t0)
        self.init_error(ode_inte.t,ode_inte.y)
        self.PID_Control(ode_inte.y[self.num_dim])
    
        time=np.array([ode_inte.t])
        x=np.array([ode_inte.y])
        u=np.array([self.u])        
        e=np.array([[self.e, self.e_dot]])
        y=np.array([self.get_y(ode_inte.y[0:self.num_dim])])
        y_ref=np.array([self.y_ref(ode_inte.t)]) 
        #print(u)
        #time.append(self.ode.t)
        #x.append(self.ode.y)
        
        while ode_inte.successful() and ode_inte.t < tf:
            ode_inte.integrate(ode_inte.t+delta_t)
            self.update_error(ode_inte.t,ode_inte.y)
            self.update_pid_control(ode_inte.y[self.num_dim])
            
            time=np.append(time,[ode_inte.t],axis=0)
            x=np.append(x,[ode_inte.y],axis=0)
            u=np.append(u,[self.u],axis=0)
            e=np.append(e,[[self.e, self.e_dot]],axis=0)
            y=np.append(y,[self.get_y(ode_inte.y[0:self.num_dim])],axis=0)
            y_ref=np.append(y_ref,[self.y_ref(ode_inte.t)],axis=0)
            #print("%g %g %g %g" % (self.ode.t, self.ode.y[0],self.ode.y[1],self.ode.y[2]))
        
        #print(e)
        self.gen_plot(time,x,u,e,y,y_ref)
        
        with open('D:\\Syno\\PID-Tuner\\ModelTesting\\'+self.name+ '.csv', 'w', newline='') as csvfile:
            fieldnames = ['t','u','y','r']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for i in range(0,len(time)):
                writer.writerow({'t':time[i],'u':u[i],'y':y[i],'r':y_ref[i]})    


        #print(plt.get_legend_handles_labels())
    def init_error(self,t,x):
        x_s=x[0:self.num_dim]
        self.e = self.y_ref(t)-self.get_y(x_s)
        self.e_dot = np.asscalar(np.array([0]))
    def update_error(self,t,x):
        x_s=x[0:self.num_dim]
        #e_i=x[self.num_dim]
        self.e = self.y_ref(t)-self.get_y(x_s)
        self.e_dot = (self.e - self.e_last)/self.dt 
        self.e_last = self.e
        #self.e_dot = (self.get_y(x_s) - self.e_last)/self.dt 
        #self.e_last = self.get_y(x_s)
        
    def update_pid_control(self,e_i):
        u = self.KP*self.e+e_i*self.KI+self.KD*self.e_dot
        self.u = np.clip(u, self.u_min, self.u_max)

    def gen_plot(self,time,x,u,e,y,y_ref):
        plt.figure(1)
        plt.subplot(221)
        plot = plt.plot(time,x)
      
        lengend_str=[ "x_%g" % (i) for i in range(1,self.num_dim+1)]
        lengend_str.append('error integral')       
        plt.legend(plot, lengend_str)
        
        plt.subplot(222)        
        plot = plt.plot(time,y,time,y_ref)     
        plt.legend(plot, ['y','y_ref'])
        
        plt.subplot(223)
        plot = plt.plot(time,u)
        plt.legend(plot,'u')
        
        plt.subplot(224)
        plot = plt.plot(time,e)
        plt.legend(plot,['e','e_dot'])
        
        plt.show()
    
class LIN_Siso_Model(Siso_Model):
    def get_y(self,x):
        return np.asscalar(self.C.dot(x))+self.u*self.D
        
    def f(self,t,x):        
        x_s=x[0:self.num_dim]
        #e_i=x[self.num_dim]  
        
        x_dot = self.A.dot(x_s) + np.transpose(self.B*self.u)       
        
        return np.concatenate((x_dot,[[self.e]]),axis=1)

     

    def __init__(self,y_ref,PID_PARA,lin_model_para):
        
        self.A=lin_model_para['A']    
        self.B=lin_model_para['B']    
        self.C=lin_model_para['C']    
        self.D=lin_model_para['D']    
        self.num_dim=self.A.shape[0]
              
        super(LIN_Siso_Model, self).__init__(y_ref,PID_PARA,lin_model_para['name'])
        
class NONLIN_Siso_Model(Siso_Model):
    def get_y(self,x):
        return np.asscalar(self.C.dot(x))+self.u*self.D
        
    def f(self,t,x,arg):        
        x_s=x[0:self.num_dim]
        #e_i=x[self.num_dim]  
        x_dot = self.A.dot(x_s) + np.transpose(B*self.u)       
        
        return np.concatenate((x_dot,[[self.e]]),axis=1)

    def __init__(self,num_dim,x_dot,y_fun,y_ref,PID_PARA):
            
        self.x_dot=x_dot
        self.y_fun=y_fun
        self.num_dim=num_dim        
        super(LIN_Siso_Model, self).__init__(y_ref,PID_PARA)        
        
        
def y_ref_t(t):
    return t         
    
def y_ref_sin(t):
    return np.sin(t)  

def y_ref_step(t):
    if t>.1:
        return 1
    else:
        return 0

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
#lin.integrate(y0=[0,0,0],t0=0,tf=10,delta_t=0.001)

lin_model_para=dict()
lin_model_para['name']  = 'model_2'
lin_model_para['A']     = np.matrix([[0,1,0],[0,0,1],[-1,-1,-1]])
lin_model_para['B']     = np.matrix([[0],[0],[1]])
lin_model_para['C']     = np.matrix([[1,1,1]])
lin_model_para['D']     = 0
    

PID_PARA=dict()
PID_PARA['P']=10
PID_PARA['I']=5
PID_PARA['D']=1
lin = LIN_Siso_Model(y_ref=y_ref_step,PID_PARA=PID_PARA,lin_model_para=lin_model_para) 
lin.integrate(y0=[1,1,1,0],t0=0,tf=5,delta_t=0.001)




















         