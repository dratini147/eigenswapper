import numpy as np
import qutip as qt

# Previous code for clarity: rho here is the solution of the Lindblad Master equation in density operator form.
# The spectral decomposition of rho is taken at all times.

# def decomposition(rho, dim):  # spectral decomposition of the density operator for each time
#     r_ev = np.zeros((tsteps, dim), dtype="float64")
#     r_es = np.empty((tsteps, dim), dtype="object")
#     for i in range(tsteps):
#         herm = rho[i].isherm
#         trace = rho[i].tr()
#         if not herm:
#             print("Density operator is not Hermitian")
#         if np.abs(trace - 1) > 0.001:
#             print("Trace is not 1")
#         ev, es = rho[i].eigenstates()
#         r_ev[i] = ev
#         r_es[i] = es
#     return r_ev, r_es


# r_ev is the array of eigenvalues (floats)
# r_es is the array of eigenvectors in QuTiP quantum object format
# The orthogonality of eigenvectors at adjacent timesteps should be close to 0.
#   When the orthogonality is high, check the next entry in the eigenvector array for a "swap".
def eigensort(r_ev, r_es, tsteps, dim):  # tracks and removes crossings in eigenvalues of the eigensolver
    for i in range(tsteps - 1):
        for j in range(dim):
            orthog = r_es[i, j].dag() * r_es[i + 1, j]
            orthog = np.abs(np.real(orthog.tr()))
            if i > 0 and np.abs(r_ev[i, j] - r_ev[i, j - 1]) < 0.1:
                orthog2 = None
                if j < dim - 1:
                    orthog2 = r_es[i, j].dag() * r_es[i + 1, j + 1]
                    orthog2 = np.abs(np.real(orthog2.tr()))
                    if orthog2 > orthog:
                        r_es[(i + 1):, [j, j + 1]] = r_es[(i + 1):, [j + 1, j]]
                        r_ev[(i + 1):, [j, j + 1]] = r_ev[(i + 1):, [j + 1, j]]
    return r_ev, r_es