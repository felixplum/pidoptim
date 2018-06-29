import numpy as np 
import casadi as ca

class PidTuner:

	def __init__(self, modelFitter, dt):
		# get order, then call .predict() with 2 args
		self.modelFitter = modelFitter

	# Objective fct: Depends on unknown y-states, known reference
	def getObjective(self, states, ref):
		obj = 0.
	def getConstraints(self, states, pid_params):
		# Consistency constraint on states
		# u_k = PID(e_k, e_k-1, itg(.), Kp, Kd, Ki)
		# y_k = f(u_k.., y_k-1..)
		print("empty")

	def fitPid(self):
		print("lol")


