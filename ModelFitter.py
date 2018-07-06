import casadi as ca
import numpy as np 
import os

class ModelFitter: 
	def __init__(self, n_input_lag, n_output_lag):
		self.n_input_lag = n_input_lag
		self.n_output_lag = n_output_lag
		self.u_coeff = None
		self.y_coeff = None
	
	def predict(self, u_prev, y_prev, u_coeff=None, y_coeff=None):
		if (u_coeff is None or y_coeff is None):
			u_coeff = self.u_coeff
			y_coeff = self.y_coeff
		assert(y_prev.shape[0] == y_coeff.shape[0])		
		assert(u_prev.shape[0] == u_coeff.shape[0])
		out = 0
		for i in range(y_prev.shape[0]):
			out += y_coeff[i]*y_prev[i]
		for i in range(u_prev.shape[0]):
			out += u_coeff[i]*(u_prev[i])
		return out

	# Set objective based on measured controls and outputs,
	# depending on optim. variables "params"
	def getObjective(self, u_in, y_out, u_coeff, y_coeff):
		obj = 0		
		max_lag = max(self.n_output_lag, self.n_input_lag)
		y_prev=[y_out[i] for i in range(0,self.n_output_lag)]
		for i in range(max_lag, len(y_out)):
			u_prev = u_in[i-self.n_input_lag:i+1] # including current input
			#y_prev = y_out[i-self.n_output_lag:i] # only past values
			y_pred = self.predict(u_prev, np.array(y_prev), u_coeff, y_coeff)
			y_prev[0:self.n_output_lag-1]=y_prev[1:self.n_output_lag]
			y_prev[self.n_output_lag-1]=y_pred
			error = y_pred - y_out[i]
			obj += error*error
		# obj += ca.dot(u_coeff, u_coeff) + 0.0001*ca.dot(y_coeff, y_coeff)
		return obj


	def fitModel(self,  u_in, y_out):
		y_coeff = ca.SX.sym('y_c', self.n_output_lag)
		u_coeff = ca.SX.sym('u_c', self.n_input_lag+1) # include current input
		params = ca.vertcat(u_coeff, y_coeff)
		obj = self.getObjective(u_in, y_out, u_coeff, y_coeff)
		# set up solver
		nlp = {'x': params, 'f': obj}
		solver_opts = {'ipopt': {'print_level': 0, 'linear_solver': 'mumps'}, 'print_time' : 0}
		solver = ca.nlpsol('solver', 'ipopt', nlp, solver_opts)
		# solver.generate_dependencies('nlp.c')
		# os.system("gcc -fPIC -shared nlp.c -o nlp.so")
		# abspath = os.path.abspath('nlp.so')
		# solver = ca.nlpsol('solver', 'ipopt', abspath, solver_opts)
		x0 = np.zeros(params.shape[0])
		res = np.array(solver(x0=x0)['x'])
		self.u_coeff = res[:u_coeff.shape[0]]
		self.y_coeff = res[u_coeff.shape[0]:]
		print(res) #
		# def eval_fct(u_prev, y_prev):
		# 	return self.predict(u_prev, y_prev, res[:u_coeff.shape[0]], res[u_coeff.shape[0]:])