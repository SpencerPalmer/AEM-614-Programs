import numpy as np
import json
import math
import matplotlib.pyplot as plt

#json read in
jsoninfo = "wingW2_json.json"

#read in json string
json_string = open(jsoninfo).read()

#Create input dictionary
input_dict = json.loads(json_string)

#Get values for everything in the dictionary
Number = input_dict["wing"]["nodes_per_semispan"]
CLa= input_dict["wing"]["airfoil_lift_slope"]
Ra = input_dict["wing"]["planform"]["aspect_ratio"]
Rt = input_dict["wing"]["planform"]["taper_ratio"]
alpha_deg = input_dict["condition"]["alpha_root[deg]"]
view = input_dict["view"]["planform"]
wing_type = input_dict["wing"]["planform"]["type"]
viewWashout = input_dict["view"]["washout_distribution"]
WashoutDistribution = input_dict["wing"]["washout"]["distribution"]
CL_design = input_dict["wing"]["washout"]["CL_design"]
TwistMagnitude = input_dict["wing"]["washout"]["magnitude[deg]"]

#Get total number of nodes N
N = Number*2-1

#Get the angle of attack, alpha in radians. Set alphaL0=0
alpha = alpha_deg*math.pi/180
alphaL0 = 0

#Get twist magnitude in radians
if (isinstance(TwistMagnitude, float)):    
    twist = TwistMagnitude*math.pi/180
    Omega = twist

#Initialize array of zeros for every array that will be eventually filled
theta = np.zeros([N])
C = np.zeros([N,N])
cTheta = np.zeros([N])
onesVector = np.ones([N])
Zb = np.zeros([N])
cb = np.zeros([N])
wTheta = np.zeros([N])

#Get the values of theta. They are equally spaced from 0-pi
for j in range(0,N):
    theta[j] = (theta[j] + j*math.pi/(N-1))

#Get c/b, use if statements to seperate if the wing is tapered or elliptic
if wing_type == "tapered":
    for i in range(0,N):
        cTheta[i] = 2/(Ra*(1+Rt))*(1-(1-Rt)*abs(math.cos(theta[i])))
elif wing_type == "elliptic":
    for i in range(0,N):
        cTheta[i] = 4/(math.pi*Ra)*math.sin(theta[i])

#Get the z/b values
for j in range(0,N):
    Zb[j]= -.5*math.cos(theta[j]) 

#Get the C matrix and an for the scenario of no twist
if WashoutDistribution == "none":
    for i in range(1,N-1):
        for j in range(0,N):
            C1 = (4/(CLa*cTheta[i]))+((j+1)/math.sin(theta[i]))
            C2 = math.sin((j+1)*theta[i])
            C[0,j]=(j+1)**2
            C[i,j] = C1*C2
            C[N-1,j]=((-1)**(j))*(j+1)**2

    #Invert the C matrix to get the an values
    an = np.linalg.inv(C).dot(onesVector)

#Get wtheta, C matrix, an, and bn for the scenario of linear twist
elif WashoutDistribution == "linear":
    for i in range(0,N):
        wTheta[i] = abs(math.cos(theta[i]))
    for i in range(1,N-1):
        for j in range(0,N):
            C1 = (4/(CLa*cTheta[i]))+((j+1)/math.sin(theta[i]))
            C2 = math.sin((j+1)*theta[i])
            C[0,j]=(j+1)**2
            C[i,j]=C1*C2
            C[N-1,j]=((-1)**(j))*(j+1)**2         
    #Invert the C matrix to get the an values
    an = np.linalg.inv(C).dot(onesVector)
    bn = np.linalg.inv(C).dot(wTheta) 

#Get wtheta, C matrix, an, and bn for the scenario of optimum twist
elif WashoutDistribution == "optimum":
    for i in range(0,N):
        if cTheta[i] < .0001:
            cTheta[i] = 1e-6
    for i in range(0,N):
        wTheta[i] = 1-(math.sin(theta[i])/(cTheta[i]/cTheta[Number-1]))
    for i in range(1,N-1):
         for j in range(0,N):
            C1 = (4/(CLa*cTheta[i]))+((j+1)/math.sin(theta[i]))
            C2 = math.sin((j+1)*theta[i])
            C[0,j]=(j+1)**2
            C[i,j]=C1*C2
            C[N-1,j]=((-1)**(j))*(j+1)**2
    #Invert the C matrix to get the an values
    an = np.linalg.inv(C).dot(onesVector)
    bn = np.linalg.inv(C).dot(wTheta) 
  
#Get KD
KD = 0
for j in range(1,N):
    KD = KD + (j+1)*(an[j]**2/an[0]**2)

#Get KL, make sure order of operations is correct
K1 = 1-(1+(math.pi*Ra/CLa))*an[0]
K2 = (1+(math.pi*Ra/CLa))*an[0]
KL = K1/K2
#Print Kappa values
print("Lift Slope Factor, kappa_L: ",KL)
print("Induced Drag Factor, kappa_D: ",KD)

#Get es
es = 1/(1+KD)
print("Span efficiency factor, e_s: ",es)

#Get the wing lift slope
CLA = CLa/((1+CLa/(math.pi*Ra))*(1+KL))
print("Wing Lift Slope, CL,a: ",CLA)

# Get CL and CDi for no twist
if WashoutDistribution == "none":
    #Get the lift coefficient
    CL = CLA*(alpha-alphaL0)
    print("Lift Coefficient, CL: ",CL)
    
    #Get the induced drag coefficient
    CDi = CL**2/(math.pi*Ra*es)
    print("Induced Drag Coefficient, CDi: ",CDi)

#Get additional Kappa values, CL, and CDi for case of linear twist
elif WashoutDistribution == "linear":
    epsilonOmega = bn[0]/an[0]
    KDL = 0
    for j in range(1,N):
        KDL = KDL+(j+1)*(an[j]/an[0])*((bn[j]/bn[0])-(an[j]/an[0]))
    KDL = KDL * 2*bn[0]/an[0]
    KDOmega = 0
    for j in range(1,N):
        KDOmega = KDOmega+(j+1)*((bn[j]/bn[0])-(an[j]/an[0]))**2
    KDOmega =  KDOmega * (bn[0]/an[0])**2
    print("epsilon_Omega: ",epsilonOmega)
    print("kappa_DL: ",KDL)
    print("kappa_DOmega: ",KDOmega)
    
    CL = CLA*((alpha - alphaL0)-epsilonOmega*Omega)
    print("Lift Coefficient, CL: ",CL)
    CDi = (CL**2*(1+KD)-KDL*CL*CLA*Omega+KDOmega*(CLA*Omega)**2)/(math.pi*Ra)
    print("Induced Drag Coefficient, CDi: ",CDi)
    
#Get additional Kappa values, CL, and CDi for case of optimum twist
elif WashoutDistribution == "optimum":
    epsilonOmega = bn[0]/an[0]
    KDL = 0
    for j in range(1,N):
        KDL = KDL+(j+1)*(an[j]/an[0])*((bn[j]/bn[0])-(an[j]/an[0]))
    KDL = KDL * 2*bn[0]/an[0]
    KDOmega = 0
    for j in range(1,N):
        KDOmega = KDOmega+(j+1)*((bn[j]/bn[0])-(an[j]/an[0]))**2
    KDOmega =  KDOmega * (bn[0]/an[0])**2
    print("epsilon_Omega: ",epsilonOmega)
    print("kappa_DL: ",KDL)
    print("kappa_DOmega: ",KDOmega)
    #Get Omega for optimum twist
    if TwistMagnitude == "optimum":
        Omega = (KDL*CL_design)/(2*KDOmega*CLA)
    CL = CLA*((alpha-alphaL0)-epsilonOmega*Omega)
    print("Lift Coefficient, CL: ",CL)
    #CDi = ((CL**2)/math.pi*Ra)*(1 + KD - ((KDL**2)/(4*KDOmega))*(2-CL_design/CL)*(CL_design/CL))
    CDi = (CL**2*(1+KD)-KDL*CL*CLA*Omega+KDOmega*(CLA*Omega)**2)/(math.pi*Ra)
    print("Induced Drag Coefficient, CDi: ",CDi)

#Plot the planform
if view == True:    
    axes = plt.gca()
    axes.plot(Zb[:],.25*cTheta[:],"r")
    axes.plot(Zb[:],-.75*cTheta[:],"r")
    for j in range(0,N):
        plt.vlines(Zb[j],-.75*cTheta[j],.25*cTheta[j], "b")
    plt.hlines(0,Zb[0],Zb[N-1], "r")
    axes.set_title("Planform")
    axes.set_xlabel("z/b")
    axes.set_ylabel("c/b")
    axes.set_aspect("equal", adjustable = "box")
    plt.show()
    
#Plot the washout distribution
if viewWashout == True:
    if WashoutDistribution == "linear":
        axes = plt.gca()
        axes.plot(Zb[:],wTheta[:],"r")
        axes.set_title("Washout Distribution")
        axes.set_xlabel("z/b")
        axes.set_ylabel("Omega")
        #axes.set_aspect("equal", adjustable = "box")
        plt.show()
    elif WashoutDistribution == "optimum":
        axes = plt.gca()
        axes.plot(Zb[:],wTheta[:],"r")
        axes.set_title("Washout Distribution")
        axes.set_xlabel("z/b")
        axes.set_ylabel("Omega")
        #axes.set_aspect("equal", adjustable = "box")
        plt.show()
    elif WashoutDistribution =="none":
        axes=plt.gca()
        plt.hlines(0,Zb[0],Zb[N-1],"r")
        axes.set_title("Washout Distribution")
        axes.set_xlabel("z/b")
        axes.set_ylabel("Omega")
        #axes.set_aspect("equal", adjustable = "box")
        plt.show()

#write the results of the C matrix, C inverse matrix, and an values to "solution.txt" txt file
if WashoutDistribution == "none":
    with open("Solution.txt", 'w') as output_file:
        output_file.write("C_matrix:\n")
        for i in range(N):
            for j in range(N):
                out_string = "{:>20.12f}".format(C[i,j])
                output_file.write(out_string)
            output_file.write("\n")
        output_file.write("C_Inverse:\n")
        for i in range(N):
            for j in range(N):
                out_string = "{:>20.12f}".format(np.linalg.inv(C)[i,j])
                output_file.write(out_string)
            output_file.write("\n")
        output_file.write("Fourier Coefficients {an}:\n")
        for i in range(N):
            out_string = "{:>20.12f}".format(an[i])
            output_file.write(out_string)
            
#Write the results for the linear washout distribution to "solution.txt"
elif WashoutDistribution == "linear":
    with open("Solution.txt", 'w') as output_file:
        output_file.write("C_matrix:\n")
        for i in range(N):
            for j in range(N):
                out_string = "{:>20.12f}".format(C[i,j])
                output_file.write(out_string)
            output_file.write("\n")
        output_file.write("C_Inverse:\n")
        for i in range(N):
            for j in range(N):
                out_string = "{:>20.12f}".format(np.linalg.inv(C)[i,j])
                output_file.write(out_string)
            output_file.write("\n")
        output_file.write("Fourier Coefficients {an}:\n")
        for i in range(N):
            out_string = "{:>20.12f}".format(an[i])
            output_file.write(out_string)
        output_file.write("\n")
        output_file.write("Fourier Coefficients {bn}:\n")
        for i in range(N):
            out_string = "{:>20.12f}".format(bn[i]) 
            output_file.write(out_string)
        output_file.write("\n")  
        
#Write the results for the optimum washout distribution to "solution.txt"
elif WashoutDistribution == "optimum":
    with open("Solution.txt", 'w') as output_file:
        output_file.write("C_matrix:\n")
        for i in range(N):
            for j in range(N):
                out_string = "{:>20.12f}".format(C[i,j])
                output_file.write(out_string)
            output_file.write("\n")
        output_file.write("C_Inverse:\n")
        for i in range(N):
            for j in range(N):
                out_string = "{:>20.12f}".format(np.linalg.inv(C)[i,j])
                output_file.write(out_string)
            output_file.write("\n")
        output_file.write("Fourier Coefficients {an}:\n")
        for i in range(N):
            out_string = "{:>20.12f}".format(an[i])
            output_file.write(out_string)
        output_file.write("\n")
        output_file.write("Fourier Coefficients {bn}:\n")
        for i in range(N):
            out_string = "{:>20.12f}".format(bn[i]) 
            output_file.write(out_string)
        output_file.write("\n")  
