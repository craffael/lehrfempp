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
err_o1 = data[:,3]
err_o2 = data[:,4]

p1_meshwidth = np.polyfit(np.log(meshwidth), np.log(err_o1), 1)
p2_meshwidth = np.polyfit(np.log(meshwidth), np.log(err_o2), 1)
p1_ndofs = np.polyfit(np.log(ndofs_o1), np.log(err_o1), 1)
p2_ndofs = np.polyfit(np.log(ndofs_o2), np.log(err_o2), 1)

plt.loglog(meshwidth, np.exp(p1_meshwidth[0]*np.log(meshwidth) + p1_meshwidth[1]), '--', label='slope: '+str(round(p1_meshwidth[0], 2)))
plt.loglog(meshwidth, np.exp(p2_meshwidth[0]*np.log(meshwidth) + p2_meshwidth[1]), '--', label='slope: '+str(round(p2_meshwidth[0], 2)))
plt.loglog(meshwidth, err_o1, 'o', label='Linear FE')
plt.loglog(meshwidth, err_o2, 'o', label='Quadratic FE')
plt.legend(loc='lower left')
plt.grid()
plt.xlabel('Mesh width [log]')
plt.ylabel('Discretization error [log]')
plt.gca().invert_xaxis()
plt.savefig(os.path.join(output_folder, 'Q_H1_meshwidth.eps'))

plt.gcf().clear()

plt.loglog(ndofs_o1, np.exp(p1_ndofs[0]*np.log(ndofs_o1) + p1_ndofs[1]), '--', label='slope: '+str(round(p1_ndofs[0], 2)))
plt.loglog(ndofs_o2, np.exp(p2_ndofs[0]*np.log(ndofs_o2) + p2_ndofs[1]), '--', label='slope: '+str(round(p2_ndofs[0], 2)))
plt.loglog(ndofs_o1, err_o1, 'o', label='Linear FE')
plt.loglog(ndofs_o2, err_o2, 'o', label='Quadratic FE')
plt.legend(loc='lower left')
plt.grid()
plt.xlabel('DOFs [log]')
plt.ylabel('Discretization error [log]')
plt.savefig(os.path.join(output_folder, 'Q_H1_dofs.eps'))
