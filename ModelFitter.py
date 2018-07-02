import casadi as ca
import numpy as np 
from csv_rw import read_model_csv
import matplotlib.pyplot as plt
import math

def get_binom_coeff_m(order):
	
	matrix = np.zeros((order+1,order+1))
	for order_cur in range(order,-1,-1):
		for row in range(order-order_cur,order+1):			
			matrix[row,order-order_cur]=(-1)**(order-row)*math.factorial(order_cur)/(math.factorial(order-row)*math.factorial(order_cur-order+row))

	
	return matrix

class ModelFitter: 
	def __init__(self, order):
		self.order =  order
		self.bin_coeff_m=get_binom_coeff_m(order)		


	
	# reutrns y_k as fct. of unknown coeff.	
	def predict(self,y_prev,u_prev,u_coeff,y_coeff):
		out = 0
		for i in range(0,self.order):
			out += u_prev[i]*u_coeff[i]
			out += -y_prev[i]*y_coeff[i]

		out = out/y_coeff[self.order]
			
			
		return out
	# Set objective based on measured controls and outputs,
	# depending on optim. variables "params"
	def getObjective(self,y_out,u_in,u_coeff,y_coeff):
		obj = 0

		y_prev=[y_out[i] for i in range(0,self.order)]
			

		
		#max_lag = max(self.n_output_lag, self.n_input_lag)
		#for i in range(max_lag, len(y_out)):
		for i in range(self.order, len(y_out)):
			u_prev=u_in[i-self.order+1:i+1]	
			

			y_pred = self.predict(y_prev,u_prev,u_coeff,y_coeff)

			error = y_pred - y_out[i]
			obj += error*error
			y_prev[0:self.order-1]=y_prev[1:self.order]
			y_prev[self.order-1]=y_pred
			
			#u_prev = u_in[i-1]
			#r_prev 		= r_in[i-3:i+1]			
			#y_prev = y_out[i-3:i] # only past values
			#y_pred = self.predict(y_prev,r_prev)
			

			
			#y0=y1
			#y1=y2
			#y2=y_pred
			#obj += error*error
			#y_prev = y_out[i-2:i+1]
			#r_prev 	= r_in[i-2:i+1]
			#e0=r_in[i-2]-y_out[i-2]
			#e1=r_in[i-1]-y_out[i-1]
			#e2=r_in[i]-y_out[i]
			#u_pred = self.get_u_2(e0,e1,e2,u0)
			#error = u_pred- u_in[i]
			#u0=u_pred
			#obj += error*error
		#obj += (alpha[0]*alpha[0]+alpha[1]*alpha[1]+beta[0]*beta[0]+beta[1]*beta[1])*0.0000001
		return obj

	def get_u_coeff(self,beta,dt):
		coeff = [0 for i in range(0,self.order)]
		
		for i in range(0,self.order):
			for p in range(0,1+i):				
				coeff[i]+=beta[self.order-1-p]*dt**(1+p)*self.bin_coeff_m[1+i,p+1]

		return coeff

	def get_u_coeff_CL(self,beta,dt,PID_coeff):
		coeff_OL = self.get_u_coeff(beta,dt)

		coeff = [0 for i in range(0,self.order+2)]
		
		for i in range(0,self.order):
			for p in range(0,3):				
				coeff[p+i]+=coeff_OL[i]*PID_coeff[p]             

		return coeff

	def get_y_coeff(self,alpha,dt):
		coeff = [0 for i in range(0,self.order+1)]
		coeff[0] =self.bin_coeff_m[0,0]	
		for i in range(1,self.order+1):
			for p in range(0,i+1):				
				if p == 0:
					coeff[i]+=self.bin_coeff_m[i,p]
				else:								
					coeff[i]+=-alpha[self.order-p]*dt**(p)*self.bin_coeff_m[i,p]	

		return coeff

	def get_y_coeff_CL(self,alpha,dt,u_coeff_CL):
		coeff_OL = self.get_y_coeff(alpha,dt)

		coeff = [0 for i in range(0,self.order+2)]		
		#coeff_pid=self.get_PID_coeff(PID_coeff,dt)
		
		for i in range(0,self.order+1):
			coeff[i]+=-coeff_OL[i]
			coeff[i+1]+=coeff_OL[i]
			coeff[i]+=u_coeff_CL[i]
		
		coeff[self.order+1]+=u_coeff_CL[self.order+1]
		

		return coeff

	

	def sim_CL(self,r,u_coeff_CL,y_coeff_CL):
		y=[0 for i in range(0,len(r))]
		#print(len(r),len(u_coeff_CL),len(y_coeff_CL))
		for i in range(len(u_coeff_CL),len(r)):
			y[i]=self.predict2(y[i-len(u_coeff_CL)+1:i],r[i-len(u_coeff_CL):i+1],u_coeff_CL,y_coeff_CL)

		return y
	def sim(self,y_out,u_in,t):
		y=[0 for i in range(0,len(u_in))]
		#print(len(r),len(u_coeff_CL),len(y_coeff_CL))

		for i in range(len(self.u_coeff),len(u_in)):
			y[i]=self.predict(y[i-len(self.u_coeff):i],u_in[i-len(self.u_coeff):i+1],self.u_coeff,self.y_coeff)

		plt.figure(1)
		plt.subplot(211)
		plot=plt.plot(t,y,t,y_out)
		plt.legend(plot,['y_fitt','y'])	
		plt.subplot(212)
		plot=plt.plot(t,y_out-y)
		plt.legend(plot,['fit_error'])	
		plt.show()
		#return y	
		

	def fitModel(self, y_out,u_in ,t):
		self.dt=np.mean(np.diff(t))

		alpha = ca.MX.sym('alpha', self.order)
		beta = ca.MX.sym('beta', self.order) # include current input
		

		params = ca.vertcat(alpha, beta)#,PID_coeff)

		u_coeff=self.get_u_coeff(beta,self.dt)
		y_coeff=self.get_y_coeff(alpha,self.dt)
		

		obj = self.getObjective(y_out,u_in,u_coeff,y_coeff)
		# set up solver
		nlp = {'x': params, 'f': obj}
		solver_opts = {'ipopt': {'print_level': 0, 'linear_solver': 'mumps'}, 'print_time' : 0}
		solver = ca.nlpsol('solver', 'ipopt', nlp, solver_opts)
		x0 = np.zeros(params.shape[0])
		res = np.array(solver(x0=x0)['x'])

		self.poles = res[0:self.order]
		self.zeros = res[self.order:]

		self.u_coeff=self.get_u_coeff(self.zeros,self.dt)
		self.y_coeff=self.get_y_coeff(self.poles,self.dt)

		print('alphas\n',self.poles)
		print('betas\n',self.zeros)
		self.sim(y_out,u_in,t)
		#def eval_fct(r):
		#	return self.sim_CL(r,u_coeff_CL,y_coeff_CL)

		#def eval_cont(u):	
		#	return self.sim_OL(u,u_coeff,y_coeff)
		#def u_func(e):	
		#	return self.sim_controller(e,self.get_PID_coeff(res2,dt))
		#return eval_fct,eval_cont,u_func

	# def __call__(self):

#in_lag= 2 
#out_lag=2
#a = ModelFitter(2) # input lag / output lag
#t,u,y,r = read_model_csv('D:\\Syno\\PID-Tuner\\ModelTesting\\model_1.csv')

#dt = np.mean(t[1:len(t)]- t[0:len(t)-1])
#print(dt)
#model,open_loop, controller = a.getFittedModel(y, r,u, dt)
#a.getFittedModel(y, r,u, dt)
#y_test=open_loop(u)
#plt.figure(1)
#plt.subplot(311)
#plot=plt.plot(t,y_test,t,y)
#plt.legend(plot,['y_fitt','y'])

#y_opt=model(r)
#plt.subplot(312)
#plot2=plt.plot(t,y_opt,t,r,t,y)
#plt.legend(plot2,['y_opt','ref','y'])

#u_out=controller(r-y_opt)
#plt.subplot(313)
#plot2=plt.plot(t,u_out)
#plt.legend(plot2,['u_out_opt'])

#y_test=np.zeros(len(t))
#u_test=np.zeros(len(t))
#y_test[0]=y[0]
#y_test[1]=y[1]
#y_test[2]=y[2]
#u_test[0]=y_test[0]
#u_test[1]=y_test[1]
#u_test[2]=y_test[2]
#for i in range(2,len(t)):
#	y_test[i]=model(y_test[i-2],y_test[i-1],u[i-1],u[i-0])
#for i in range(3,len(t)):	
#	u_test[i] =controller(u_test[i-3],u_test[i-2],u_test[i-1],r[i-3],r[i-2],r[i-1],r[i]) 

#y_test =[ model(np.array([y[t_c-3],y[t_c-2],y[t_c-1]]),np.array([r[t_c-3],r[t_c-2],r[t_c-1],r[t_c-0]])) for t_c in range(3,len(t))]

 

#u_test =[controller(np.array([r[t_c-2],r[t_c-1],r[t_c-0]])-np.array([y[t_c-2],y[t_c-1],y[t_c-0]]),u[t_c-1]) for t_c in range(2,len(t))]
#model = a.getFittedModel(y, u)
#y_test =[ model(np.array([y[t_c-1], y[t_c]]), np.array([u[t_c-3],u[t_c-2],u[t_c-1]])) for t_c in range(out_lag,len(t))]




#
#
#

#plt.show()
#print("Prediction is: ", model(np.array([u[2], u[3]]), np.array([y[0],y[1],y[2]])))