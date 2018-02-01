#! /usr/bin/env python3
# -*- coding: utf8 -*-

##########################################################################
#                                                                        #
# Programa: irc4orca                        Date: 29/01/2018             #
#                                                                        #
# Usage: irc4orca.py inpfilename                                         #
#                                                                        #
# This is a wrapper that implements Morokuma's IRC algorithm in Cartesian#
# coordinates (J. Phys. Chem. 66, 2153), with some modifications.        #
#                                                                        #
# The input file is a normal Orca input file, although I recommend it not#
# to bear any indication as to the type of calculation (Opt, enGrad, SP  #
# or MD keywords). The parameters to be read by the wrapper are given as #
# comments, one instruction per line:                                    #
#                                                                        #
# Mandatory instructions:                                                #
#                                                                        #
# #orcacmd                    - location of an executable file that      #
#                               takes an orca input as argument and      #
#                               outputs Orca's input as a filen with     #
#                               the same base name and the.out extension #
#                                                                        #
# Mandatory instructions (unless it's a restart):                        #
#                                                                        #
# #irchess filename.hess      - location of the file hess file for the   #
#                               system                                   #
#                                                                        #
# #ircmode n                  - Vibrational mode to follow (default:0)   #
#                                                                        #
# Optional instructions:                                                 #
#                                                                        #
# #ircrestart [0/1]           - This is a restart from a point in the    #
#                               IRC that is not the TS (default:0/False) #
#                                                                        #
# #ircguess filename.gbw      - location of a GBW for the first point    #
#                                                                        #
# #ircalpha x.xx              - Alpha parameter for scaling the initial  #
#                               displacement (default: 0.1).             #
#                                                                        #
# #ircdelta x.xx              - Scaling parameter for finding the closest#
#                               minimum at a given point (default: 0.05).#
#                                                                        #
# #ircdamp x.xx               - Percentage of the previous displacement  #
#                               on the calculation of the current one    #
#                               (default: 0.05).                         #
#                                                                        #
# #ircautodamp [0/1]          - Update ircdamp. The update procedure     #
#                               reduces ircdump by 0.1/log(n) for n>1    #
#                               (default: 0, False)                      #
#                                                                        #
# #irctol x.xx                - Acceptable value of the RMS gradient for #
#                               considering the IRC to have converged    #
#                               (default:1.0e-4).                        #
#                                                                        #
# #ircdir [-1/+1]             - Direction of the initial displacement    #
#                               (default: +1).                           #
#                                                                        #
# #ircmaxd                    - Maximum norm of the displacement vector  #
#                               (default: 0.01 angs/amu^(1/2)).          #
#                                                                        #
# #ircpts                     - Maximum number of points in the IRC      #
#                               (default: 25)                            #
#                                                                        #
# #ircalg [1/2]               - IRC algorithm: 1 (default) is the        #
#                               traditional Morokuma algorithm.          #
#                               2 is an updated version that requires    #
#                               3 additional energy evaluations per step.#
#                                                                        #
#   NOTES: ircalpha and ircmaxd may be used together to fine tune the    #
#          development of the IRC procedure. Small values of ircmaxd     #
#          tend to make the calculation stop near the TS, so a larger    #
#          ircalpha may be desired. Another solution to the problem of   #
#          starting an IRC with small ircmaxd is to increase ircdamp     #
#          up to 0.9 and use ircautodamp to reduce it when the gradient  #
#          is string enough to carry the geometry towards its minimum.   #
#                                                                        #
#          The calculation will stop when any of the following is        #
#          verified: ircpts is reached; the gradient is bellow irctol,   #
#          or the energy of the current geometry in the IRC is higher    #
#          than that of the previous point.                              #
#                                                                        #
#          The procedure updates the ircdamp parameter by multiplying by #
#          0.1/log(n), in which n is the cycle number. This starts for   #
#          n=2 and ircdamp is automatically set to zero when it goes     #
#          bellow 1.0e-5.                                                #
#                                                                        #
#                                                                        #
# (c) Filipe Teixeira 2012 (Current version:08/12/2013)                  #
#                                                                        #
##########################################################################

# No editing is required bellow this line #

import os
import sys
import numpy as np

#####################################
# Defining classes to store data in #
#####################################

class Atom():
	def __init__(self, line):
		l=line.split()
		self.symbol=l[0]
		self.coords=np.array(list(map(float,l[1:])))
	def printxyz(self):
		return "%s %9.6f %9.6f %9.6f\n"%(self.symbol, self.coords[0], self.coords[1], self.coords[2])

class ToolKit():
	# just a placeholder for all common data
	def __init__(self,name):
		natoms=1
		self.alpha=0.1
		self.basename='.'.join(name.split('.')[:-1])
		self.out=open("%s.log"%(self.basename),'a',1)
		self.delta=0.05
		self.direction=1
		self.energies=[]
		self.energy=0.00
		self.geos=[]
		self.restart=False
		self.geometry=[]
		self.guessfn=''
		self.damp=0.05
		self.algorithm=1
		self.autodamp=False
		self.prevgrad=0.0
		self.hessfn=''
		self.maxdispl=0.01
		self.mode=0
		self.dgrad=0.0
		self.natoms=natoms
		self.npoints=25
		self.template=[]
		self.tolerance=1.0e-04
		self.orcacmd='UNDEFINED'
		self.ReadInput(name)
		self.ReadHessian()
	def printPars(self):
		stmp="  %13s: %s\n"
		ftmp="  %13s: %7.4f\n"
		itmp="  %13s: %7d\n"
		self.out.write('')
		self.out.write(itmp%('Algorithm',self.algorithm))
		self.out.write(itmp%('N. Points',self.npoints))
		self.out.write(ftmp%('Grad. Tol.',self.tolerance))
		self.out.write('')
		self.out.write(stmp%('Hessian',self.hessfn))
		self.out.write(itmp%('Mode',self.mode))
		self.out.write(itmp%('Direction',self.direction))
		self.out.write('')
		self.out.write(ftmp%('Alpha',self.alpha))
		self.out.write(ftmp%('Delta',self.delta))
		self.out.write('')
		self.out.write(ftmp%('Damp Factor',self.damp))
		self.out.write(itmp%('Damp Update',self.autodamp))
		self.out.write('')
		self.out.write(ftmp%('Max. Displ.',self.maxdispl))
		self.out.write('')
		self.out.write(stmp%('Guess',self.guessfn))
		self.out.write("\n------------------------------------------------\n")
	def ReadInput(self,name):
		#parse inp file
		ifile=open(name,'r')
		idata=ifile.readlines()
		ifile.close()
		ingeo=False
		for line in idata:
			copy=True
			if ingeo:
				copy=False
			if '#' in line:
				copy=False
				if 'irches' in line.lower(): #hess file name
					l=line.split()
					self.hessfn=l[-1]
				elif 'orcacmd' in line.lower(): # orca executable
					l=line.split()
					self.orcacmd=l[-1]
				elif 'ircguess' in line.lower(): # gbw file name for the initial guess
					l=line.split()
					self.guessfn=l[-1]
				elif 'ircalpha' in line.lower(): # alpha
					l=line.split()
					self.alpha=float(l[-1])
				elif 'ircdelta' in line.lower(): # delta
					l=line.split()
					self.delta=float(l[-1])
				elif 'ircdamp' in line.lower(): #damping x(prevdisp)+(1-x)disp
					l=line.split()
					self.damp=float(l[-1])
				elif 'ircrestart' in line.lower(): #update damping?
					l=line.split()
					if (int(l[-1])==1):
						self.restart=True
					else:
						self.restart=False
				elif 'ircautodamp' in line.lower(): #update damping?
					l=line.split()
					if (int(l[-1])==1):
						self.autodamp=True
					else:
						self.autodamp=False
				elif 'irctol' in line.lower(): 
					l=line.split()
					self.tolerance=float(l[-1])
				elif 'ircalg' in line.lower(): # algorithm
					l=line.split()
					self.algorithm=int(l[-1])
				elif 'ircdir' in line.lower(): # direction +1 or -1
					l=line.split()
					self.direction=int(l[-1])
				elif 'ircmode' in line.lower(): # mode number (same numbering as in hess)
					l=line.split()
					self.mode=int(l[-1])
				elif 'ircmaxd' in line.lower(): # number of points to do
					l=line.split()
					self.maxdispl=float(l[-1])
				elif 'ircpts' in line.lower(): # number of points to do
					l=line.split()
					self.npoints=int(l[-1])
			if ('*' in line):
				if (ingeo):
					ingeo=False
					copy=False
				else:
					ingeo=True
			if (('*' not in line) and ingeo):
				self.geometry.append(Atom(line))
			if (copy):
				line=line.lower()
				line=line.replace('sp',' ')
				line=line.replace('engrad',' ')
				line=line.replace('opt',' ')
				line=line.replace('numfreq',' ')
				line=line.replace('moread',' ')
				if ('moinp' in line): 
					line='\n'
				self.template.append(line)
		self.natoms=len(self.geometry)
	def ReadHessian(self):
		#read hessian from hess file 
		hfile=open(self.hessfn,'r')
		hdata=hfile.readlines()
		hfile.close()
		modes=np.zeros((3*self.natoms,3*self.natoms))
		energy=0.00
		modestart=0
		massstart=0
		for i in range(len(hdata)):
			if '$act_energy' in hdata[i]:
				self.energies.append(float(hdata[i+1]))
			if '$normal_modes' in hdata[i]:
				modestart=i+2
			if '$atoms' in hdata[i]:
				massstart=i+2
		done = False
		i=modestart
		cols=[0,1,2,3,4,5]
		while (not done):
			if '#' in hdata[i]:
				done=True
			elif (hdata[i]=='\n'):
				done=True
			elif ('.' in hdata[i]):
				l=hdata[i].split()
				line=int(l[0])
				l=list(map(float,l[1:]))
				modes[line,cols]=l
			elif (hdata[i].strip != ''):
				cols=list(map(int,hdata[i].split()))
			i=i+1
		done = False
		i=massstart
		j=0
		self.mass=np.zeros(self.natoms)
		while (not done):
			if '#' in hdata[i]:
				done=True
			elif (hdata[i]=='\n'):
				done=True
			elif ('.' in hdata[i]):
				l=hdata[i].split()
				self.mass[j]=float(l[1])
				j += 1
			i=i+1
		self.displacement=modes[:,self.mode] 
		self.grad=np.zeros(3*self.natoms) #allocating space for the gradients

#######################
# Interface with Orca #
#######################

	
def doEnergy(geo,pars):
	#run energy calculation -> returns energy
	name=pars.basename+'.tmp.sp'
	#os.system("rm %s*"%name) #clean previous calculations with the same name
	inpfile=open(name+'.inp','w')
	if (pars.guessfn==""):
		inpfile.write("! SP\n")
	else:
		inpfile.write("! SP MoRead\n%%moinp \"%s\"\n\n"%(pars.guessfn))
	for line in pars.template:
		inpfile.write(line)
	for a in geo:
		inpfile.write(a.printxyz())
	inpfile.write("*\n\n")
	inpfile.close()
	os.system("%s %s.inp > %s.out"%(pars.orcacmd, name, name))
	opipe=os.popen("grep 'FINAL SINGLE POINT ENERGY' %s.out"%(name))
	odata=opipe.readlines()
	opipe.close()
	os.system("rm %s*"%name) #clean up
	return float(odata[-1].split()[-1])

def doGrad(geom,pars):
	#run gradient calculation -> returns energy and maxgrad
	# also updates pars.grad vector
	#os.system("rm %s*"%(name)) #clean previous calculations with the same name
	name=pars.basename+'.tmp.grd'
	inpfile=open(name+'.inp','w')
	newguess=pars.basename+".i4o.last.gbw"
	if (pars.guessfn==""):
		inpfile.write("! EnGrad\n")
	else:
		inpfile.write("! EnGrad MoRead\n%%moinp \"%s\"\n\n"%(pars.guessfn))
	for line in pars.template:
		inpfile.write(line)
	for a in geom:
		inpfile.write(a.printxyz())
	inpfile.write("*\n\n")
	inpfile.close()
	os.system("%s %s.inp"%(pars.orcacmd, name))
	os.system("mv %s.gbw %s"%(name,newguess))
	pars.guessfn=newguess
	ofile=open("%s.out"%(name),'r')
	odata=ofile.readlines()
	ofile.close()
	ingrad=False
	atom=0
	for line in odata:
		if (ingrad) and ':' in line:
			ldata=line.split()
			ldata=list(map(float,ldata[-3:]))
			s=3*atom
			e=s+3
			pars.grad[s:e]=ldata
			atom=atom+1
		if 'CARTESIAN GRADIENT' in line:
			ingrad=True
		if (atom>=pars.natoms):
			break
	op=os.popen("grep 'RMS gradient' %s.out"%(name))
	od=op.readlines()
	op.close()
	maxgrad=float(od[-1].split()[-1])
	op=os.popen("grep 'FINAL SINGLE POINT ENERGY' %s.out"%(name))
	od=op.readlines()
	op.close()
	energy=float(od[-1].split()[-1])
	os.system("rm %s*"%(name)) #clean previous calculations with the same name
	return (energy, maxgrad)

######################################
# Geometry manipulation and printing #
######################################

def geodisplace(geo, pars, dvec):
	#adjust geo for displacement dvec
	displace=dvec
	newgeo=[]
	for i in geo:
		newgeo.append(Atom(i.printxyz()))
	i=0
	for a in newgeo:
		a.coords[0] += displace[i]
		a.coords[1] += displace[i+1]
		a.coords[2] += displace[i+2]
		i += 3
	return newgeo

def printTrj(params,n):
	trj=open(params.basename+'.trj','a')
	trj.write("%d\nIRC for Orca point %d E=%14.7f\n"%(params.natoms,n,params.energy))
	for a in params.geometry:
		trj.write(a.printxyz())
	trj.close()

##############
# IRC kernel #
##############

def Morokuma(pars,start=False):
	#Does a cycle in the Morokuma algorithm, returns the energy and MaxGrad of the
	# new point
	if (start):
		#apply the appropriate sign to the displacement and convert to angs.
		pars.displacement=float(pars.direction)*pars.displacement
		#for i in range(len(pars.mass)):
		#	pars.displacement[(3*i)] /= np.sqrt(pars.mass[i])
		#	pars.displacement[(3*i)+1] /= np.sqrt(pars.mass[i])
		#	pars.displacement[(3*i)+2] /= np.sqrt(pars.mass[i])
		pars.displacement=pars.alpha*(pars.displacement/np.linalg.norm(pars.displacement))
		#displace geometry following the normal mode
		geo1=geodisplace(pars.geometry,pars,pars.displacement)
		E1,MG1=doGrad(geo1,pars)
		#generate vector D, adapting from eq 6 in J. Chem. Phys. 66, 2153
		tmpvec=pars.grad
		D=-(pars.displacement/np.linalg.norm(pars.displacement))+(tmpvec/np.linalg.norm(tmpvec))
		#D=(pars.displacement/np.linalg.norm(pars.displacement))-(tmpvec/np.linalg.norm(tmpvec))
	else:
		#scale down gradients and calculate the displacement
		pars.displacement=(pars.damp*pars.displacement)-((1.0-pars.damp)*pars.grad)
		for i in range(len(pars.mass)):
			pars.displacement[(3*i)] /= np.sqrt(pars.mass[i])
			pars.displacement[(3*i)+1] /= np.sqrt(pars.mass[i])
			pars.displacement[(3*i)+2] /= np.sqrt(pars.mass[i])
		pars.displacement=pars.maxdispl*(pars.displacement/np.linalg.norm(pars.displacement))
		#displace geometry following the gradient
		geo1=geodisplace(pars.geometry,pars,pars.displacement)
		E1,MG1=doGrad(geo1,pars)
		#generate vector D, adapting from eq 6 in J. Chem. Phys. 66, 2153
		tmpvec=pars.grad
		D=(pars.displacement/np.linalg.norm(pars.displacement))-(tmpvec/np.linalg.norm(tmpvec))
	# find the optimum delta that will assure the new point is at the 
	# local minimum
	if (pars.algorithm==1):
		geo2=geodisplace(geo1,pars,pars.delta*D)
		E2=doEnergy(geo2,pars)
		if (E2>E1):
			newdelta=0.5*pars.delta
		else:
			newdelta=2.0*pars.delta
		geo3=geodisplace(geo1,pars,newdelta*D)
		E3=doEnergy(geo3,pars)
		Evals=np.array([E1, E2, E3])
		Deltavals=np.array([0.0,pars.delta,newdelta])
	else:
		geo2=geodisplace(geo1,pars,pars.delta*D)
		E2=doEnergy(geo2,pars)
		delta3=0.5*pars.delta
		geo3=geodisplace(geo1,pars,delta3*D)
		E3=doEnergy(geo3,pars)
		delta4=2.0*pars.delta
		geo4=geodisplace(geo1,pars,delta4*D)
		E4=doEnergy(geo4,pars)
		delta5=-0.5*pars.delta
		geo5=geodisplace(geo1,pars,delta5*D)
		E5=doEnergy(geo5,pars)
		Evals=np.array([E1, E2, E3, E4, E5])
		Deltavals=np.array([0.0,pars.delta,delta3, delta4, delta5])
	deltaFit=np.polyfit(Deltavals,Evals,deg=2)
	opdelta=-(deltaFit[1]/(2.0*deltaFit[0]))
	# update the geometry to the new point, update and return energy and gradients
	pars.geos.append(pars.geometry)
	pars.energies.append(pars.energy)
	pars.geometry=geodisplace(geo1,pars,opdelta*D)
	newE, newMaxGrad = doGrad(pars.geometry,pars) 
	pars.energy=newE
	#the following line is just for debugging purposes
	#print(Evals, Deltavals,np.linalg.norm(pars.displacement), opdelta)
	return (newE, newMaxGrad)

def ircdrv(inpname):
	params=ToolKit(inpname)
	params.out.write("""   IRC wrapper for Orca - version 2.0
   by Filipe Teixeira, 
   REQUIMTE
   Faculdade de Ciencias da Universidade do Porto

  Using Orca from: %s 
  
  Parameters for this run: \n"""%(params.orcacmd))
	params.printPars()
	keep=True
	n=0
	oldE=0.0
	params.out.write("   Pt. %20s %9s %10s\n"%('Energy','RMS Grad.','Damp'))
	params.out.write("------------------------------------------------\n")
	while (keep):
		if (params.restart or (n>0)):
			E,MG=Morokuma(params,False)
		else:
			E,MG=Morokuma(params,True)
		n=n+1
		if(params.autodamp and (n>1)):
			params.damp = params.damp * (0.1/np.log(n))
			if (params.damp<1.0e-5):
				params.damp=0.0
				params.autodamp=False
		params.out.write('@  %3d %20.9f %9.5f %10.2e\n'%(n, E, MG, params.damp))
		if (n>params.npoints):
			keep=False
			printTrj(params,n) #appends newGeo
			params.out.write("----------------------------------------------\n")
			params.out.write("--          IRC CALCULATION ENDED           --\n")
			params.out.write("--             MAXPTS ACHIEVED!             --\n")
			params.out.write("----------------------------------------------\n")
		elif (MG<params.tolerance):
			keep=False
			printTrj(params,n) #appends newGeo
			params.out.write("----------------------------------------------\n")
			params.out.write("--          IRC CALCULATION ENDED           --\n")
			params.out.write("--              TOL ACHIEVED!               --\n")
			params.out.write("----------------------------------------------\n")
		elif (E>oldE):
			keep=False
			#printTrj(params,n) #appends newGeo
			params.out.write("----------------------------------------------\n")
			params.out.write("--             ENERGY INCREASED             --\n")
			params.out.write("--            GEOMETRY  MIGHT BE            --\n")
			params.out.write("--          VERY CLOSE TO A MINIMUM         --\n")
			params.out.write("--                                          --\n")
			params.out.write("--        IRC CALCULATION TERMINATED!       --\n")
			params.out.write("----------------------------------------------\n")
		else:
			printTrj(params,n) #appends newGeo
			oldE=E
	params.out.close()



if __name__=='__main__':
	if(len(sys.argv)!=2):
		print("""IRC4Orca - Version 2.0
An Implementation of Morokuma's IRC method for the Orca ESS Software.
by Filipe Teixeira

Usage: %s file.inp

Please consult https://github.com/teixeirafilipe/irc4orca for more information.

"""%(sys.argv[0]))
		sys.exit(1)
	ircdrv(sys.argv[1])

