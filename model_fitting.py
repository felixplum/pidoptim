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
		
	def get_u(self,e,u):
   		out = 0
   		
   		out = self.h1*e[0] + self.h2*e[1] + self.h3*e[2] + u
   		return out

	def get_u_2(self,e,u,PID_coeff):
		out = 0   		
		out += PID_coeff[0]*e[0] + PID_coeff[1]*e[1] + PID_coeff[2]*e[2] + u
		return out
	# reutrns y_k as fct. of unknown coeff.
	def predict(self,y_prev,r_prev):
		v1=y_prev[0]-y_prev[1]
		v2=y_prev[1]-y_prev[2]
		v3=y_prev[2]

		out = 0
		out += ( v1*self.gamma1+v2*self.gamma2+v3*self.gamma3 + self.p1*(self.h1*(-y_prev[1]+r_prev[1])+ self.h2*(-y_prev[2]+r_prev[2])+self.h3*r_prev[3])+\
			self.p2*(self.h1*(v1+r_prev[1]-r_prev[0])+self.h2*(v2+r_prev[2]-r_prev[1])+self.h3*(v3+r_prev[3]-r_prev[2])))/self.gamma0
			
		return out

	def predict2(self,y_prev,u_prev,u_coeff,y_coeff):
		out = 0
		
		for i in range(0,self.order+1):
			out += u_prev[i]*u_coeff[i]
			out += -y_prev[i]*y_coeff[i]
		out += u_prev[self.order+1]*u_coeff[self.order+1]	
		out = out/y_coeff[self.order+1]
			
			
		return out

	def predict3(self,y_prev,u_prev,u_coeff,y_coeff):
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
			

			y_pred = self.predict3(y_prev,u_prev,u_coeff,y_coeff)

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

	def getObjective_PID(self,y_out,ref,u_coeff_CL,y_coeff_CL,PID_coeff):
		obj = 0
		g = 0
		y_prev=[y_out[i] for i in range(0,self.order+1)]
			
		u0=0
		
		#max_lag = max(self.n_output_lag, self.n_input_lag)
		#for i in range(max_lag, len(y_out)):
		for i in range(self.order+1, len(y_out)):
			
			r_prev=ref[i-self.order-1:i+1]	
			

			y_pred = self.predict2(y_prev,r_prev,u_coeff_CL,y_coeff_CL)

			error = ref[i]- y_pred
			obj += error*error*1000
			y_prev[0:self.order]=y_prev[1:self.order+1]
			y_prev[self.order]=y_pred

		
			
			#y_prev = y_out[i-3:i] # only past values
			#y_pred = self.predict2(y0,y1,y2,r_prev)
			#error = r_in[i]- y_pred
			#y0=y1
			#y1=y2
			#y2=y_pred
			#obj += error*error

			e=r_prev[1:]-y_prev			
			u_pred = self.get_u_2(e,u0,PID_coeff)
			u0=u_pred			
			error = u_pred
			obj += error*error*0.01
			#u1=u_pred
			g =ca.vertcat(g,u_pred)

		#obj += 0.0001*(PID_coeff[0]*PID_coeff[0]+PID_coeff[1]*PID_coeff[1]+PID_coeff[2]*PID_coeff[2])
		return obj,g

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

	def get_PID_coeff(self,PID_coeff,dt):
		coeff = [0 for i in range(0,3)]
		coeff[0]=(PID_coeff[2])/dt
		coeff[1]=(-dt*PID_coeff[0]+PID_coeff[1]*dt*dt*0.5-2*PID_coeff[2])/dt
		coeff[2]=(dt*dt*PID_coeff[1]*0.5+PID_coeff[0]*dt+PID_coeff[2])/dt
		return coeff

	def sim_CL(self,r,u_coeff_CL,y_coeff_CL):
		y=[0 for i in range(0,len(r))]
		#print(len(r),len(u_coeff_CL),len(y_coeff_CL))
		for i in range(len(u_coeff_CL),len(r)):
			y[i]=self.predict2(y[i-len(u_coeff_CL)+1:i],r[i-len(u_coeff_CL):i+1],u_coeff_CL,y_coeff_CL)

		return y
	def sim_OL(self,u,u_coeff,y_coeff):
		y=[0 for i in range(0,len(u))]
		#print(len(r),len(u_coeff_CL),len(y_coeff_CL))

		for i in range(len(u_coeff),len(u)):
			y[i]=self.predict3(y[i-len(u_coeff):i],u[i-len(u_coeff):i+1],u_coeff,y_coeff)

		return y	
	def sim_controller(self,e,PID_coeff):
		y=[0 for i in range(0,len(e))]
		#print(len(r),len(u_coeff_CL),len(y_coeff_CL))

		for i in range(3,len(e)):
			y[i]=self.get_u_2(e[i-2:i+1],y[i-1],PID_coeff) 

		return y	

	def getFittedModel(self, y_out, ref,u_in , dt):
		alpha = ca.MX.sym('alpha', self.order)
		beta = ca.MX.sym('beta', self.order) # include current input
		PID_coeff = ca.MX.sym('pid_c', 3) 

		params = ca.vertcat(alpha, beta)#,PID_coeff)


		u_coeff=self.get_u_coeff(beta,dt)
		y_coeff=self.get_y_coeff(alpha,dt)




		#self.dt=dt


		

		#self.gamma1=1
		#self.gamma2=-2+dt*alpha[1]
		#self.gamma3=1-dt*dt*alpha[0]-dt*alpha[1]
		#self.gamma0=2-dt*dt*alpha[0]-dt*alpha[1]+(dt*beta[0]+beta[1])*self.h3

		#self.p1=(beta[0]*dt+beta[1])*dt
		#self.p2=beta[1]*dt

		

		

		obj = self.getObjective(y_out,u_in,u_coeff,y_coeff)
		# set up solver
		nlp = {'x': params, 'f': obj}
		solver_opts = {'ipopt': {'print_level': 0, 'linear_solver': 'mumps'}, 'print_time' : 0}
		solver = ca.nlpsol('solver', 'ipopt', nlp, solver_opts)
		x0 = np.zeros(params.shape[0])
		res = np.array(solver(x0=x0)['x'])
		#print(res) #


		u_coeff_CL=self.get_u_coeff_CL(res[self.order:],dt,self.get_PID_coeff(PID_coeff,dt))
		y_coeff_CL=self.get_y_coeff_CL(res[0:self.order],dt,u_coeff_CL)
		
		

		obj,g = self.getObjective_PID(y_out,ref,u_coeff_CL,y_coeff_CL,self.get_PID_coeff(PID_coeff,dt))
		
		params = ca.vertcat(PID_coeff)
		
		nlp = {'x': params, 'f': obj,'g':g}
		solver_opts = {'ipopt': {'print_level': 0, 'linear_solver': 'mumps'}, 'print_time' : 0}
		solver = ca.nlpsol('solver', 'ipopt', nlp, solver_opts)
		x0 = np.zeros(params.shape[0])
		res2 = np.array(solver(x0=x0,lbg=-400,ubg=400)['x'])
		
		u_coeff=self.get_u_coeff(res[self.order:],dt)
		y_coeff=self.get_y_coeff(res[0:self.order],dt)

		u_coeff_CL=self.get_u_coeff_CL(res[self.order:],dt,self.get_PID_coeff(res2,dt))
		y_coeff_CL=self.get_y_coeff_CL(res[0:self.order],dt,u_coeff_CL)

		print('poles\n',res[0:self.order])
		print('zeros\n',res[self.order:])
		print('K_p ',res2[0][0])
		print('K_i ',res2[1][0])
		print('K_d ',res2[2][0])
		def eval_fct(r):
			return self.sim_CL(r,u_coeff_CL,y_coeff_CL)

		def eval_cont(u):	
			return self.sim_OL(u,u_coeff,y_coeff)
		def u_func(e):	
			return self.sim_controller(e,self.get_PID_coeff(res2,dt))
		return eval_fct,eval_cont,u_func

	# def __call__(self):

#in_lag= 2 
#out_lag=2
a = ModelFitter(4) # input lag / output lag
t,u,y,r = read_model_csv('D:\\Syno\\PID-Tuner\\ModelTesting\\model_2.csv')

dt = np.mean(t[1:len(t)]- t[0:len(t)-1])
#print(dt)
model,open_loop, controller = a.getFittedModel(y, r,u, dt)
y_test=open_loop(u)
plt.figure(1)
plt.subplot(311)
plot=plt.plot(t,y_test,t,y)
plt.legend(plot,['y_fitt','y'])

y_opt=model(r)
plt.subplot(312)
plot2=plt.plot(t,y_opt,t,r,t,y)
plt.legend(plot2,['y_opt','ref','y'])

u_out=controller(r-y_opt)
plt.subplot(313)
plot2=plt.plot(t,u_out)
plt.legend(plot2,['u_out_opt'])

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


plt.show()
#print("Prediction is: ", model(np.array([u[2], u[3]]), np.array([y[0],y[1],y[2]])))