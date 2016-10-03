import math
class Metrics:
	def __init__(self,t,p):
		self.RMSE = 0.0
		self.AAE = 0.0
		self.PEARSON = 0.0
		if(len(t)!=len(p)):
			exit("TARGET AND PREDICTED LIST HAVE DIFFERENT LENGTH\nlen(t)"+str(len(t))+" len(p)"+str(len(p)))
		self.t=t
		self.p=p
		self.N=len(t)
	def computeRMSE(self):
		temp = 0.0
		for i in range(0,self.N):
			temp += (self.t[i]-self.p[i])*(self.t[i]-self.p[i])
		temp = math.sqrt(temp/float(self.N))
		self.RMSE = temp
	def computeAAE(self):
		temp = 0.0
		for i in range(0,self.N):
			temp += abs(self.t[i]-self.p[i])
		temp = temp/float(self.N)
		self.AAE = temp
	def computeStdAAE(self):
		stdAAE = 0.0
		for i in range(0,self.N):
			stdAAE+=math.pow((abs(self.t[i]-self.p[i])-self.AAE),2)
		stdAAE = math.sqrt(stdAAE/float((self.N)-1))
		self.stdAAE = stdAAE
		
	def computeAverage(self):
		avT = 0.0
		avP = 0.0
		for i in range(0,self.N):
			avT+=self.t[i]
			avP+=self.p[i]
		avT = avT/float(self.N)
		avP = avP/float(self.N)
		self.avT = avT
		self.avP = avP

	def computeStd(self):
		self.computeAverage()		
		stdT = 0.0
		stdP = 0.0
		for i in range(0,self.N):
			stdT+=math.pow((self.t[i]-self.avT),2)
			stdP+=math.pow((self.p[i]-self.avP),2)
		stdT = math.sqrt(stdT/float((self.N)-1))
		stdP = math.sqrt(stdP/float((self.N)-1))
		self.stdT = stdT
		self.stdP = stdP
						
	def computePEARSON(self):
		r = 0.0
		#average t
		aT = 0.0
		#average p
		aP = 0.0
		#st
		st = 0.0
		#sp
		sp = 0.0
		#calculating average t and p
		for i in range(0,self.N):
			aT += self.t[i]
			aP += self.p[i]
		aT = aT/float(self.N)
		aP = aP/float(self.N)
		############################
		#calculating st and sp
		for i in range(0,self.N):
			st += (self.t[i]-aT)*(self.t[i]-aT)
			sp += (self.p[i]-aP)*(self.p[i]-aP)
		
		######################
		#calculating r
		for i in range(0,self.N):
			r += self.t[i]*self.p[i]
		
		r = (r-float(self.N)*aT*aP)/math.sqrt(st*sp)
		r = r*r
		self.PEARSON = r	
		

