import numpy as np
import matplotlib.pyplot as plt

threshold_iter = 1000

print("Enter a value for length of the channel L:")
L = input()
print("Enter a value for height of the channel H:")
H = input()
print("Enter a value for h1")
h1 = input()
print("Enter a value for h2")
h2 = input()
print("Enter a value for N1")
N1 = input()
print("Enter a value for N2")
N2 = input()
print("Enter a value for Vin")
vin = input()

# L = 15.0
# H = 10.0
# h1 = 2.0
# h2 = 4.0
# N1 = 301
# N2 = 201

# vin=4	


vout = (float)(H*vin)/(h2-h1)
delX = 	float(L)/(N1)
delY =  float(H)/(N2)

beta_val = ((L/float(N1-1) **2)/(H/float(N2-1)**2))
beta_divide = 2*(1+beta_val)


def calculate_velocities(N1, N2, phi):
	velocity_xcomp = np.zeros((N2, N1))  # To store all the x values of the velocity
	velocity_ycomp = np.zeros((N2, N1))  # To store all the y values of the velocity
	for i in range(N2):
	    for j in range(N1):
	        if j - 1 >= 0:
	            velocity_xcomp[i][j] = (phi[i,j] - phi[i,j - 1]) / delX
	        else:
	            velocity_xcomp[i][j] = vin
	        if i - 1 >= 0:
	            velocity_ycomp[i][j] = (phi[i,j] - phi[i-1,j]) / delY
	        else:
	            velocity_ycomp[i][j] = 0
	return velocity_xcomp, velocity_ycomp


# Function to set all the boundary conditions after every iteration
def setBoundaryConditions(phi):

	# for IN surface
	for i in range(0,N2):
		phi[i,0]=1

	# for TOP surface
	for j in range(0,N1):
		phi[N2-1,j]=phi[N2-2,j]

	# for BOTTOM surface
	for j in range(0,N1):
		phi[0,j]=phi[1,j]

	#for RIGHT1 surface
	for i in range(0,int(N2 * (1 - h1 / H))):
		phi[i,N1-1]=phi[i,N1-2]

	#for RIGHT2 surface
	for i in range(int(N2 * (1 - h2 / H)),N2):
		phi[i,N1-1]=phi[i,N1-2]

	# for OUT surface
	for i in range(int(N2 * (1 - h2 / H)) , int(N2 * (1 - h1 / H))):
		phi[i,N1-1]=phi[i,N1-2] + delX*vout

	return phi



phi = np.random.rand(N2, N1)   #Initialised phi matrix with random values
for i in range(0,N2):
	for j in range(0, N1):
		phi[i, j] = j*delX


setBoundaryConditions(phi)     # set the boundaries conditions initially


# Iterative solution starts
print("Iterations start")
for iteration in range(0, threshold_iter):
	print("iteration: %d" % iteration)
	for i in range(1, N2-1):
		for j in range(1, N1-1): 
			phi[i,j] =  (phi[i,j-1] + phi[i,j+1] + beta_val*phi[i+1][j] + beta_val*phi[i-1][j])/float(beta_divide)
	phi = setBoundaryConditions(phi)
			#print("returned back")

print("Iterations Finished")
# print(phi)


velocity_xcomp, velocity_ycomp = calculate_velocities(N1, N2, phi)


# Saving the phis and the velocities to a file.
np.savetxt('/home/sanidhya/Desktop/Semester_6/Continuum Mechanics/phi.txt', phi, fmt='%.18e', delimiter=' ', newline='\n')
np.savetxt('/home/sanidhya/Desktop/Semester_6/Continuum Mechanics/Xvelocity.txt', velocity_xcomp, fmt='%.18e', delimiter=' ', newline='\n')
np.savetxt('/home/sanidhya/Desktop/Semester_6/Continuum Mechanics/Yvelocity.txt', velocity_ycomp, fmt='%.18e', delimiter=' ', newline='\n')


# to plot the contour plots
x, y = np.arange(0, L, delX), np.arange(0, H, delY)
x, y = np.meshgrid(x, y)
fig = plt.figure()
colorinterpolation = 50
colorMap = plt.cm.spectral
plt.contourf(x, y, phi, colorinterpolation, cmap=colorMap)
plt.title('Contour Plot for phis')
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.show()


# To plot the Quiver Plot
plt.figure()
plt.title('Quiver Plot for the velocity')
Q = plt.quiver(x[::4, ::4], y[::4, ::4], velocity_xcomp[::4, ::4], velocity_ycomp[::4, ::4], units='width')
plt.show()
