import numpy as np 
import casadi as ca
import matplotlib.pyplot as plt
class PidTuner:

	def __init__(self, modelFitter):
		# get order, then call .predict() with 2 args
		self.modelFitter = modelFitter
		self.order=modelFitter.order

	# Objective fct: Depends on unknown y-states, known reference
	def get_PID_coeff(self,PID_coeff,dt):
		coeff = [0 for i in range(0,3)]
		coeff[0]=(PID_coeff[2])/dt
		coeff[1]=(-dt*PID_coeff[0]+PID_coeff[1]*dt*dt*0.5-2*PID_coeff[2])/dt
		coeff[2]=(dt*dt*PID_coeff[1]*0.5+PID_coeff[0]*dt+PID_coeff[2])/dt
		return coeff

	def predict(self,y_prev,u_prev,u_coeff,y_coeff):
		out = 0
		
		for i in range(0,self.order+1):
			out += u_prev[i]*u_coeff[i]
			out += -y_prev[i]*y_coeff[i]
		out += u_prev[self.order+1]*u_coeff[self.order+1]	
		out = out/y_coeff[self.order+1]
			
			
		return out	

	def get_u(self,e,u,PID_coeff):
		out = 0   		
		out += PID_coeff[0]*e[0] + PID_coeff[1]*e[1] + PID_coeff[2]*e[2] + u
		return out

	def sim(self,y_in,r,t,PID_coeff):
		y=[0 for i in range(0,len(r))]
		#print(len(r),len(u_coeff_CL),len(y_coeff_CL))

		for i in range(len(self.u_coeff_CL),len(r)):
			y[i]=self.predict(y[i-len(self.u_coeff_CL)+1:i],r[i-len(self.u_coeff_CL)+1:i+1],self.u_coeff_CL,self.y_coeff_CL)

		
		u=self.sim_controller(r-y,PID_coeff)
		plt.figure(1)
		plt.subplot(211)
		plot=plt.plot(t,y,t,y_in,t,r)
		plt.legend(plot,['y_opt','y','r'])	
		plt.subplot(212)
		plot=plt.plot(t,u)
		plt.legend(plot,['u_opt'])	
		plt.show()	

	def sim_controller(self,e,PID_coeff):
		y=[0 for i in range(0,len(e))]
		#print(len(r),len(u_coeff_CL),len(y_coeff_CL))

		for i in range(3,len(e)):
			y[i]=self.get_u(e[i-2:i+1],y[i-1],PID_coeff) 

		return y	
	def getObjective(self,y_out,ref,u_coeff_CL,y_coeff_CL,PID_coeff):
		obj = 0
		g = 0
		y_prev=[y_out[i] for i in range(0,self.order+1)]			
		u0=0

		for i in range(self.order+1, len(y_out)):
			
			r_prev=ref[i-self.order-1:i+1]
			y_pred = self.predict(y_prev,r_prev,u_coeff_CL,y_coeff_CL)

			error = ref[i]- y_pred
			obj += error*error*1000
			y_prev[0:self.order]=y_prev[1:self.order+1]
			y_prev[self.order]=y_pred	

			e=r_prev[1:]-y_prev			
			u_pred = self.get_u(e,u0,PID_coeff)
			u0=u_pred			
			error = u_pred
			obj += error*error*0.001
		
			g =ca.vertcat(g,u_pred)

		return obj,g
	def getConstraints(self, states, pid_params):
		# Consistency constraint on states
		# u_k = PID(e_k, e_k-1, itg(.), Kp, Kd, Ki)
		# y_k = f(u_k.., y_k-1..)
		print("empty")

	def fitPid(self, y_out, ref,u_in,u_bounds,t):
		dt = self.modelFitter.dt

		PID_coeff = ca.MX.sym('pid_c', 3) 
		alpha=self.modelFitter.poles
		zeros=self.modelFitter.zeros
		u_coeff_CL=self.modelFitter.get_u_coeff_CL(zeros,dt,self.get_PID_coeff(PID_coeff,dt))
		y_coeff_CL=self.modelFitter.get_y_coeff_CL(alpha,dt,u_coeff_CL)
		
		

		#obj,g = self.getObjective(y_out,ref,u_coeff_CL,y_coeff_CL,self.get_PID_coeff(PID_coeff,dt))
		obj,g = self.getObjective(y_out,ref,u_coeff_CL,y_coeff_CL,self.get_PID_coeff(PID_coeff,dt))
		params = ca.vertcat(PID_coeff)
		
		#nlp = {'x': params, 'f': obj}
		nlp = {'x': params, 'f': obj,'g':g}
		solver_opts = {'ipopt': {'print_level': 0, 'linear_solver': 'mumps'}, 'print_time' : 0}
		solver = ca.nlpsol('solver', 'ipopt', nlp, solver_opts)
		x0 = np.zeros(params.shape[0])
		res2 = np.array(solver(x0=x0,lbg=u_bounds[0],ubg=u_bounds[1])['x'])

		self.u_coeff_CL=self.modelFitter.get_u_coeff_CL(zeros,dt,self.get_PID_coeff(res2,dt))
		self.y_coeff_CL=self.modelFitter.get_y_coeff_CL(alpha,dt,self.u_coeff_CL)



		#res2 = np.array(solver(x0=x0)['x'])

		#u_coeff=self.get_u_coeff(res[self.order:],dt)
		#y_coeff=self.get_y_coeff(res[0:self.order],dt)

		#u_coeff_CL=self.get_u_coeff_CL(res[self.order:],dt,self.get_PID_coeff(res2,dt))
		#y_coeff_CL=self.get_y_coeff_CL(res[0:self.order],dt,u_coeff_CL)

		print('K_p ',res2[0][0])
		print('K_i ',res2[1][0])
		print('K_d ',res2[2][0])

		self.sim(y_out,ref,t,self.get_PID_coeff(res2,dt))