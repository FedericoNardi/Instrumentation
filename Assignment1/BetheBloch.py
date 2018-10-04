import numpy as np
from math import pi
import matplotlib.pyplot as plt

def BetheBloch(Beam, material):
	# The function takes as input two lists, corresponding to incoming Beam and absorber material, specified as
	# material ={"A": \massNumber,"Z": \atomicNumber,"I":\MeanExcitationEnergy[MeV],"rho":\Density[g/cm3]} and
	# Beam = {"Energy":[MeV], "Particlemass":[Mev/c2], "z":, "Type": 0 for heavy particles, 1 for electrons}
	# Constants
	electronmass = 0.511 #Mev
	classicalradius = 2.817e-13 #cm

	K = 2*pi*6.022e23*classicalradius**2*electronmass

	#Relativistic parameters
	gamma = Beam["Energy"]/Beam["Particlemass"] + 1
	BetaSquare = 1 - 1/gamma**2

	if Beam["Type"]==1:
		Tmax = Beam["Energy"]
		F = 1-BetaSquare+(Tmax**2/8 - (2*Tmax + 1)*np.log(2) )/((Tmax + 1)**2)
	else:
		Tmax = 2*electronmass*BetaSquare*gamma**2 / ( 1 + 2*electronmass*gamma/Beam["Particlemass"] + (electronmass/Beam["Particlemass"])**2 )
		F = 0

	Prefactor = material["rho"]*K*material["Z"]/material["A"]*Beam["z"]**2/BetaSquare
	argument = 2*electronmass*BetaSquare*gamma**2*Tmax / (material["I"]**2)

	StoppingPower = Prefactor * ( np.log(argument) - 2*BetaSquare + F )

	return StoppingPower

#PROBLEM 1 ----------------------------------------------------------------------------------------------------------

Silicon = {"Z":14, "A":28.085, "I":173e-6, "rho":2.33} 			
Germanium = {"Z":32, "A":72.63, "I":350e-6, "rho":5.323}

AlphaBeam = {"Energy": 5.00 , "Particlemass":3.727e3, "z":2, "Type":0} 
				
print("SILICON ABSORBER - Alpha beam")
print("Energy Loss: ", BetheBloch(AlphaBeam, Silicon))#/Silicon["rho"])
print("---------------------")
print("GERMANIUM ABSORBER - Alpha beam")
print("Energy Loss: ", BetheBloch(AlphaBeam, Germanium))#/Germanium["rho"])
print("---------------------")

ElectronBeam = {"Energy": 5.00 , "Particlemass":0.511, "z":1, "Type":1}

print("SILICON ABSORBER - Electron beam")
print("Energy Loss: ", BetheBloch(ElectronBeam, Silicon))
print("---------------------")
print("GERMANIUM ABSORBER - Electron beam")
print("Energy Loss: ", BetheBloch(ElectronBeam, Germanium))
print("---------------------")



#PROBLEM 2 ----------------------------------------------------------------------------------------------------------

ProtonBeam = {"Energy": 0 , "Particlemass":938.3, "z":1, "Type":0}
Portland = {"Z":14.2, "A":28.87, "I":132.2e-6, "rho":2.300}

#Components of Portland Concrete
Hydrogen = {"Z":1,"A":1.008,"I":19.3e-6,"rho":8.37e-5,"Percentage":0.01}
Carbon = {"Z":6,"A":12.011,"I":81e-6,"rho":2.00,"Percentage":0.001}
Oxygen = {"Z":8,"A":16.00,"I":95e-6,"rho":1.33e-3,"Percentage":0.529}
Sodium = {"Z":11,"A":22.99,"I":149e-6,"rho":9.71e-1,"Percentage":0.016}
Magnesium = {"Z":12,"A":24.31,"I":156e-6,"rho":1.74,"Percentage":0.002}
Aluminum = {"Z":13., "A":26.981, "I":166e-6, "rho":2.70, "Percentage":0.034}
Silicon = {"Z":14, "A":28.085, "I":173e-6, "rho":2.33,"Percentage":0.337}
Potassium = {"Z":19, "A":39.10,"I":190e-6,"rho":8.2e-1,"Percentage":0.013}
Calcium = {"Z":20,"A":40.08,"I":191e-6,"rho":1.55,"Percentage":0.044}
Iron = {"Z":26,"A":55.84,"I":286e-6,"rho":7.87,"Percentage":0.014}

points = 500
Stoppingpower_PC_Bragg = np.zeros(points)
Stoppingpower_PC = np.zeros(points)

for i in range(points):
	E = np.linspace(0.1e3, 10e3, points)

	ProtonBeam["Energy"]=E[i]
	Stoppingpower_PC[i] = BetheBloch(ProtonBeam,Portland)
	Stoppingpower_PC_Bragg[i] = (BetheBloch(ProtonBeam,Hydrogen)*Hydrogen["Percentage"]+
		BetheBloch(ProtonBeam,Carbon)*Carbon["Percentage"]+
		BetheBloch(ProtonBeam,Oxygen)*Oxygen["Percentage"]+
		BetheBloch(ProtonBeam,Sodium)*Sodium["Percentage"]+
		BetheBloch(ProtonBeam,Magnesium)*Magnesium["Percentage"]+
		BetheBloch(ProtonBeam,Aluminum)*Aluminum["Percentage"]+
		BetheBloch(ProtonBeam,Silicon)*Silicon["Percentage"]+
		BetheBloch(ProtonBeam,Potassium)*Potassium["Percentage"]+
		BetheBloch(ProtonBeam,Calcium)*Calcium["Percentage"]+
		BetheBloch(ProtonBeam,Iron)*Iron["Percentage"])

plt.semilogx(E,Stoppingpower_PC_Bragg)
plt.semilogx(E,Stoppingpower_PC)
plt.legend(["Portland Concrete w/ Bragg rule","Portland Concrete"])
plt.title('Stopping Power: Bragg rule vs Empyrical data')
plt.xlabel("Beam Energy [MeV]")
plt.ylabel("Stopping power [MeV/cm]")
plt.grid(True, which="both", linestyle='--')
plt.savefig("figures/Bragg")
plt.show()


Stoppingpower_Al = np.zeros(points)
E = np.logspace(2, 4, points)

for i in range(points):

	ProtonBeam["Energy"]=E[i]
	Stoppingpower_PC[i] = BetheBloch(ProtonBeam,Portland)
	Stoppingpower_Al[i] = BetheBloch(ProtonBeam,Aluminum)

plt.plot(E,Stoppingpower_Al)#, marker='x',markersize=3)
plt.plot(E,Stoppingpower_PC)#, marker='x')
plt.legend(["Aluminum","Portland Concrete"])
plt.title('Stopping Power: Aluminum vs Portland Concerte')
plt.xlabel("Beam Energy [MeV]")
plt.ylabel("Stopping power [MeV/cm]")
plt.grid(True, which="both", linestyle='--')
axs=plt.gca()
axs.set_xscale('log')
plt.ylim(0,1e3)
plt.xlim(1e2,1e4)
plt.autoscale()
plt.savefig("figures/AlvsPC")
plt.show()

def BetheBloch2(Beam, material):
	# Constants
	electronmass = 0.511 #Mev
	classicalradius = 2.817e-13 #cm

	K = 2*pi*6.022e23*classicalradius**2*electronmass

	#Relativistic parameters
	gamma = Beam["Energy"]/Beam["Particlemass"] + 1
	BetaSquare = 1 - 1/gamma**2

	C = -4.24
	a = 0.0802
	m = 3.63
	X1 = 3.01
	X0 = 0.1708
	X = np.log10(np.sqrt(BetaSquare)*gamma)

	if X<X0:
		delta = 0
	else:
		if X<X1:
			delta = 4.6052*X + C + a*(X1-X)**m
		else:
			delta = 4.6052*X + C

	if Beam["Type"]==1:
		Tmax = Beam["Energy"]
		F = 1-BetaSquare+(Tmax**2/8 - (2*Tmax + 1)*log(2) )/((Tmax + 1)**2)
	else:
		Tmax = 2*electronmass*BetaSquare*gamma**2 / ( 1 + 2*electronmass*gamma/Beam["Particlemass"] + (electronmass/Beam["Particlemass"])**2 )
		F = 0

	Prefactor = material["rho"]*K*material["Z"]/material["A"]*Beam["z"]**2/BetaSquare
	argument = 2*electronmass*BetaSquare*gamma**2*Tmax / (material["I"]**2)

	StoppingPower = Prefactor * ( np.log(argument) - 2*BetaSquare + F - delta)#/material["rho"]

	return StoppingPower

points = 5000

StoppingPower = np.zeros(points)
StoppingPower2 = np.zeros(points)

E = np.logspace(2, 4, points)

for i in range(points):
	ProtonBeam["Energy"] = E[i]
	StoppingPower[i] = BetheBloch(ProtonBeam,Aluminum)
	StoppingPower2[i] = BetheBloch2(ProtonBeam,Aluminum)

plt.plot(E,StoppingPower,'r-')#, marker='x',markersize=3)
plt.plot(E,StoppingPower2,'b-')#, marker='x')
plt.legend(["Without correction","With correction"])
plt.title('Stopping power of Al - Density correction')
plt.xlabel("Beam Energy [MeV]")
plt.ylabel("Stopping power [MeV/cm]")
plt.grid(True, which="both", linestyle='--')
axs=plt.gca()
axs.set_xscale('log')
axs.set_yscale('log')


#Find minimum ionization energies
StoppingPower = np.asarray(StoppingPower)
StoppingPower2 = np.asarray(StoppingPower2)
index = StoppingPower.argmin()
index2 = StoppingPower2.argmin()
print("-------- Minimum ionization energies --------")
print("Al without correction: ",E[index],"MeV")
print("Al with correction: ",E[index2],"MeV")
print("---------------------------------------------")

plt.axvline(x=E[index],color='r',linewidth=0.5)
plt.axvline(x=E[index2],color='b',linewidth=0.5)

plt.autoscale()
plt.savefig("figures/BetheVsBethe++")
plt.show()

from scipy import integrate

Air = {"Z":0.76*7+0.29*8+0.01*18,"A":0.76*14.007+0.29*16.00+0.01*39.95,"I":85.7e-6,"rho":1.20e-3}

IntegrationPoints = 50
RangePoints = 100

InverseStoppingPower = np.zeros(IntegrationPoints)
Range = np.zeros(RangePoints)
Emax = np.linspace(0.1e3, 1e4, RangePoints)

for k in range(RangePoints):

	E = np.linspace(0, Emax[k], IntegrationPoints)

	for i in range(IntegrationPoints):
		ProtonBeam["Energy"] = E[i]
		InverseStoppingPower[i] = 1/BetheBloch(ProtonBeam,Aluminum)
	Range[k] = integrate.simps(InverseStoppingPower,E)

print("Range for a 10 MeV Beam in Aluminum: ", Range[0])

plt.loglog(Emax, Range)
plt.title('Proton beam in Aluminum - range')
plt.xlabel("Beam Energy [MeV]")
plt.ylabel("Range [cm]")
plt.grid(True, which="both", linestyle='--')
axs=plt.gca()
axs.set_xscale('log')
axs.set_yscale('log')
plt.autoscale()
plt.savefig("figures/Aluminumrange")
plt.show()

#PROBLEM 5 ----------------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
def compton(energy, angle):
	DeltaInverseEnergy = 1/0.511*(1-np.cos(angle))
	InverseFinalEnergy = 1/energy + DeltaInverseEnergy
	return 1/InverseFinalEnergy

# E = 2 #MeV
# E1 = compton(E,np.pi/6)
# E2 = compton(E,np.pi/3)

# print("Initial energy: ",E)
# print("Energy after 1st collision: ",E1)
# print("Energy after 2nd collision: ",E2)
# print("Total energy loss: ", E-E2)

# Compton edges

E0 = [0.06, 1.332, 2.164]
# points = 1000
# Angles = np.linspace(0, 2*np.pi, points)
# DepositedEnergy=np.zeros(points)
# DepositedEnergy1=np.zeros(points)
# DepositedEnergy2=np.zeros(points)

E = np.zeros(3)
Edges = np.zeros(3)
for k in range(3):
	E[k] = compton(E0[k],np.pi)
	Edges[k] = E0[k]-E[k]

for k in range(3):
	print("-----------------")
	print("Energy: ", E0[k])
	print("Backward scattered energy: ", E[k])
	print("Deposition: ", Edges[k])



# plt.plot(Angles, DepositedEnergy)

# for k in range(points):
# 	DepositedEnergy1[k] = E0[1] - compton(E0[1],Angles[k])

# plt.plot(Angles, DepositedEnergy1)

# for k in range(points):
# 	DepositedEnergy2[k] = E0[2] - compton(E0[2],Angles[k])

# plt.semilogy(Angles, DepositedEnergy2)

# plt.axvline(x=np.pi,color='r',linewidth=0.5)
# plt.axvline(x=np.pi,color='b',linewidth=0.5)
# plt.title("Compton scattering - energy deposition")
# plt.legend(["0.06 MeV","1.332 MeV","2.164 MeV"])
# plt.xlabel("\gamma angle [rad]")
# plt.ylabel("Energy [MeV]")
# plt.grid()
# plt.show()
# plt.savefig("figures/Edges")


