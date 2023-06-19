# Unprojected Augmented Krylov Methods

A library which tests and compares unprojected augmented Krylov subspace methods with the corresponding projected augmented methods
and the standard (non augmented) method. 

This code is associated to the preprint 

[1] Burke, Liam, and Kirk M. Soodhalter. "Augmented unprojected Krylov subspace methods" arXiv preprint arXiv:2206.12315 (2022).

This library has the two folders.

## Main
Contains the main files used to implement all the revelant functions:

* `fom.m`: A function which solves a linear system Ax = b using the standard restarted Full Orthogonalization Method (FOM).
* `r_fom.m`: A function which solves a linear system Ax = b using the (projected) recycled Full Orthogonalization Method (rFOM).
* `ur_fom.m`: A function which solves a linear system Ax = b using the (un projected) recycled Full Orthogonalization Method (urFOM).

* `gmres.m`: A function which solves a linear system Ax = b using the standard restarted Generalized Minimum Residual Method (GMRES).
* `r_gmres.m`: A function which solves a linear system Ax = b using the (projected) recycled Generalized Minimum Residual Method (rGMRES).
* `ur_gmres.m`: A function which solves a linear system Ax = b using the (un projected) recycled Generalized Minimum Residual Method (urGMRES).

* `gmres.m`: A function which solves a linear system Ax = b using the standard restarted Generalized Minimum Residual Method (GMRES).
* `r_gmres.m`: A function which solves a linear system Ax = b using the (projected) recycled Generalized Minimum Residual Method (rGMRES).
* `ur_gmres.m`: A function which solves a linear system Ax = b using the (un projected) recycled Generalized Minimum Residual Method (urGMRES).

* `r_fom_ritz_recycling.m`: A function which constructs a recycling subspace for rFOM using an augmented Krylov subspace K + U, by solving
  the appropriate Ritz problem for rFOM.

* `r_fom_ritz_recycling.m`: A function which constructs a recycling subspace for rFOM using an augmented Krylov subspace K + U, by solving
   the appropriate Ritz problem for rFOM.

* `ur_fom_ritz_recycling.m`: A function which constructs a recycling subspace for urFOM using an augmented Krylov subspace K + U, by solving
   the appropriate Ritz problem for urFOM.

* `ur_gmres_harmritz_recycling.m`: A function which constructs a recycling subspace for urGMRES using an augmented Krylov subspace K + U, by solving
  the appropriate harmonic Ritz problem for urGMRES.

* `r_gmres_harmritz_recycling.m`: A function which constructs a recycling subspace for rGMRES using an augmented Krylov subspace K + U, by solving
  the appropriate harmonic Ritz problem for rGMRES.

## Tests

Contains all the files used to run each of the tests

* `fom_qcd_test.m`: Tests and compares unprojected recycled FOM (ur_fom) to recycled FOM (r_fom) and standard FOM (fom) by solving a sequence of linear systems with a QCD matrix. Convergence curves and number of vectors the coefficient matrix A is applied to is recorded to compare each algorithm.

* `fom_neumann_test.m`: Tests and compares unprojected recycled FOM (ur_fom) to recycled FOM (r_fom) and standard FOM (fom) by solving a sequence of linear systems with a Neumann matrix. Convergence curves and number of vectors the coefficient matrix A is applied to is recorded to compare each algorithm.

* `gmres_qcd_test.m`: Tests and compares unprojected recycled GMRES (ur_gmres) to recycled GMRES (r_gmres) and standard GMRES (gmres) by solving a sequence of linear systems with a QCD matrix. Convergence curves and number of vectors the coefficient matrix A is applied to is recorded to compare each algorithm.

* `fom_residual_approx_test.m` and  `gmres_residual_approx_test.m` : Respectively demonstrate the use of the residual estimate norm for unprojected methods proposed in [1], by comparing the convergence curve of the residual norm estimate to the norm of the exact residual.


