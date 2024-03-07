import numpy as np
import json
import math

#json read in
jsonfile = "json_input.json"

#read in json string
json_string = open(jsonfile).read()

# create input dictionary
input_dict = json.loads(json_string)
print(input_dict)

filename = input_dict["geometry"]

with open(filename, "r") as f:
    
    #initilize info list
    info = []
    
    for line in f:
        
        #split the line
        info.append([float(point) for point in line.split() ])
        
    #close file
    f.close()

#turn into numpy array
info = np.array(info)
n = len(info)

alphaDeg = input_dict["alpha[deg]"]  
alpha = alphaDeg*np.pi/180
FreestreamVel = input_dict["freestream_velocity"]

#Arrays of zeros for values
A = np.zeros([n,n])
b = np.zeros([n])
Lj =np.zeros([n-1])
xj = np.zeros([n])
yj = np.zeros([n])
cx = np.zeros([n-1])
cy = np.zeros([n-1])
solve = np.zeros([n-1,n-1,2,2])

#Get x,y of each node
xj = info[:,0]
yj = info[:,1]

#Get x,y of control points
for i in range(0,n-1):
    cx[i] = (xj[i] + xj[i+1])/2
    cy[i] = (yj[i] + yj[i+1])/2

#Get length Lj, 1.6.19
for j in range(0, n-1):
    Lj[j] = math.sqrt((xj[j+1] - xj[j])**2 + (yj[j+1] - yj[j])**2)


for i in range(0, n-1):
    for j in range(0, n-1):
        delx = xj[j+1] - xj[j]
        dely = yj[j+1] - yj[j]
        
        #Matrix M1, 1.6.20
        M1 = np.array([[delx, dely],[-dely, delx]])
        
        #Matrix M2, 1.6.20
        M2 = np.array([[cx[i]-xj[j]],[cy[i]-yj[j]]])
        
        #Xi, Eta equation 1.6.20
        XN = (1/Lj[j])*M1.dot(M2)
        
        #Make it easier to call Xi and Eta
        X = XN[0]
        N = XN[1]
        
        #Get Psi, equation 1.6.22
        psi = .5*np.log((X**2+N**2)/((X-Lj[j])**2+N**2))
        
        #Get Phi, equation 1.6.21
        phi = np.arctan2(N*Lj[j], N**2 + X**2 - X*Lj[j])
        
        #Matrix N1, 1.6.23
        N1 = np.array([[delx, -dely],[dely, delx]])
        
        #Matrix N2, 1.6.23
        N2 = np.array([[(Lj[j] - X)*phi + N*psi, X*phi - N*psi], [N*phi - (Lj[j] -X)*psi-Lj[j],-N*phi-X*psi+Lj[j]]])
        N2 = N2.reshape(2,2)
        
        #Matrix P, 1.6.23        
        P = (1/(2*math.pi*Lj[j]**2)*N1.dot(N2))
        solve[j,i,0:2,0:2] = P
        
        #Get matrix A, Equation 1.6.25
        A[i,j] = A[i,j] + (xj[i+1]-xj[i])/Lj[i]*P[1,0]-(yj[i+1]-yj[i])/Lj[i]*P[0,0]
        A[i,j+1] = A[i,j+1] + (xj[i+1]-xj[i])/Lj[i]*P[1,1]-(yj[i+1]-yj[i])/Lj[i]*P[0,1]

#Kutta condition 1.6.26 and 1.6.27
A[n-1,0]=1
A[n-1,n-1]=1

#Get b Vector, Equation 1.6.28
for j in range(0,n-1):
    b[j] = FreestreamVel*((yj[j+1]-yj[j])*np.cos(alpha)-(xj[j+1]-xj[j])*np.sin(alpha))/Lj[j]
b[n-1]=0    

#Get gamma Vector, Equation 1.6.28, Use inverse of A matrix.
gamma = np.linalg.inv(A).dot(b)

#Get lift coefficient, 1.6.32
CL = 0
for j in range(0,n-1):
    CL = CL+Lj[j]*(gamma[j]+gamma[j+1])/FreestreamVel

#Get moment about the leading edge, 1.6.33
CMle = 0
for j in range(0,n-1):
    Xle = 2*xj[j]*gamma[j]+xj[j]*gamma[j+1]+xj[j+1]*gamma[j]+2*xj[j+1]*gamma[j+1]
    Yle = 2*yj[j]*gamma[j]+yj[j]*gamma[j+1]+yj[j+1]*gamma[j]+2*yj[j+1]*gamma[j+1]
    CMle = CMle + Lj[j]*(Xle*np.cos(alpha)/FreestreamVel+Yle*np.sin(alpha)/FreestreamVel)
CMle = CMle*(-1/3)

#Get momenet at the quarter chord, test 1 equation sheet.
CMc4 = CMle+CL/4

print([CL,CMle,CMc4])
import matplotlib.pyplot as plt
axes=plt.gca()
axes.plot(info[:,0],info[:,1])
plt.show()
#print(info[0])
#print(info[1])