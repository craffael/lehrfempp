import sys
import os
import numpy as np
import matplotlib.pyplot as plt

input_file = os.path.abspath(sys.argv[1])
output_folder = os.path.dirname(input_file)

data = np.genfromtxt(input_file, delimiter=',')
meshwidth = data[:,0]
ndofs_o1 = data[:,1]
ndofs_o2 = data[:,2]
H1_err_o1 = data[:,3]
H1_err_o2 = data[:,4]

H1_p1_meshwidth = np.polyfit(np.log(meshwidth), np.log(H1_err_o1), 1)
H1_p2_meshwidth = np.polyfit(np.log(meshwidth), np.log(H1_err_o2), 1)
H1_p1_ndofs = np.polyfit(np.log(ndofs_o1), np.log(H1_err_o1), 1)
H1_p2_ndofs = np.polyfit(np.log(ndofs_o2), np.log(H1_err_o2), 1)

plt.loglog(meshwidth, np.exp(H1_p1_meshwidth[0]*np.log(meshwidth) + H1_p1_meshwidth[1]), '--', label='slope: '+str(round(H1_p1_meshwidth[0], 2)))
plt.loglog(meshwidth, np.exp(H1_p2_meshwidth[0]*np.log(meshwidth) + H1_p2_meshwidth[1]), '--', label='slope: '+str(round(H1_p2_meshwidth[0], 2)))
plt.loglog(meshwidth, H1_err_o1, 'o', label='Linear FE')
plt.loglog(meshwidth, H1_err_o2, 'o', label='Quadratic FE')
plt.legend(loc='lower left')
plt.grid()
plt.xlabel('Mesh width [log]')
plt.ylabel('Discretization error [log]')
plt.title('Discretization errors with respect to $H^1$ semi-norm')
plt.gca().invert_xaxis()
plt.savefig(os.path.join(output_folder, 'L_H1_meshwidth.eps'))

plt.gcf().clear()

plt.loglog(ndofs_o1, np.exp(H1_p1_ndofs[0]*np.log(ndofs_o1) + H1_p1_ndofs[1]), '--', label='slope: '+str(round(H1_p1_ndofs[0], 2)))
plt.loglog(ndofs_o2, np.exp(H1_p2_ndofs[0]*np.log(ndofs_o2) + H1_p2_ndofs[1]), '--', label='slope: '+str(round(H1_p2_ndofs[0], 2)))
plt.loglog(ndofs_o1, H1_err_o1, 'o', label='Linear FE')
plt.loglog(ndofs_o2, H1_err_o2, 'o', label='Quadratic FE')
plt.legend(loc='lower left')
plt.grid()
plt.xlabel('DOFs [log]')
plt.ylabel('Discretization error [log]')
plt.title('Discretization errors with respect to $H^1$ semi-norm')
plt.savefig(os.path.join(output_folder, 'L_H1_dofs.eps'))
