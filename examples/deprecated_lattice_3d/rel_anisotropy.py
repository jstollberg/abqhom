import numpy as np

# load stiffness file
file = r"C:\Users\jonat\Documents\abqhom\examples\lattice_3d\results\result_0.01_210.0_0.3.csv"
C = np.genfromtxt(file, delimiter=",")

# compute best isotropic approximation
Q = np.array([[1, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0],
              [0, 0, 1, 0, 0, 0],
              [0, 0, 0, 2, 0, 0],
              [0, 0, 0, 0, 2, 0],
              [0, 0, 0, 0, 0, 2]], dtype=float)
A1 = 1/3*np.array([[1, 1, 1, 0, 0, 0],
                   [1, 1, 1, 0, 0, 0],
                   [1, 1, 1, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0]], dtype=float)
A2 = 1/(6*np.sqrt(5))*np.array([[4, -2, -2, 0, 0, 0],
                                [-2, 4, -2, 0, 0, 0],
                                [-2, -2, 4, 0, 0, 0],
                                [0, 0, 0, 3, 0, 0],
                                [0, 0, 0, 0, 3, 0],
                                [0, 0, 0, 0, 0, 3]], dtype=float)
C0 = (np.trace(np.dot(np.dot(np.dot(Q, C), Q), A1))*A1 
      + np.trace(np.dot(np.dot(np.dot(Q, C), Q), A2))*A2)
Cr = C - C0

# transform basis
C[0:3,3:6] *= np.sqrt(2)
C[3:6,0:3] *= np.sqrt(2)
C[3:6,3:6] *= 2
Cr[0:3,3:6] *= np.sqrt(2)
Cr[3:6,0:3] *= np.sqrt(2)
Cr[3:6,3:6] *= 2

# compute relative anisotropy
delta = np.linalg.norm(Cr)/np.linalg.norm(C)
print(delta)