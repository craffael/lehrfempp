import sys
import os
import numpy as np
import matplotlib.pyplot as plt

input_file = os.path.abspath(sys.argv[1])
output_folder = os.path.dirname(input_file)

data = np.genfromtxt(input_file, delimiter=',')
p = data[:,0]
ndofs_Q = data[:,1]
ndofs_L = data[:,2]
H1_err_Q = data[:,3]
L2_err_Q = data[:,4]
H1_err_L = data[:,5]
L2_err_L = data[:,6]

plt.semilogy(p**2, L2_err_Q, 'o', label='$L^2$ norm')
plt.semilogy(p**2, H1_err_Q, 'o', label='$H^1 semi-norm$')
plt.grid()
plt.xlabel('$p^2$')
plt.ylabel('Discretization error [log]')
plt.title('Discretization errors for p-FEM')
plt.savefig(os.path.join(output_folder, 'Square_degPoly.eps'))

plt.gcf().clear()

plt.semilogy(ndofs_Q, L2_err_Q, 'o', label='$L^2$ norm')
plt.semilogy(ndofs_Q, H1_err_Q, 'o', label='$H^1 semi-norm$')
plt.grid()
plt.xlabel('N')
plt.ylabel('Discretization error [log]')
plt.title('Discretization errors for p-FEM')
plt.savefig(os.path.join(output_folder, 'Square_N.eps'))
