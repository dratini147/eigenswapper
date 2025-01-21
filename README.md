This code was written for QuTiP (Quantum Toolbox in Python).
Matrix diagonalizers often sort eigenvalues in ascending order (I mean, how else to sort them logically?). Sometimes this is a problem, if you are tracking the eigenvalues and eigenvectors with respect to a changing parameter --
e.g. a density operator evolving over time.

This takes input arrays of eigenvalues/vectors and 'swaps' them to track continuous eigenvectors, by comparing the orthogonality of eigenvectors at adjacent time(parameter) steps.
