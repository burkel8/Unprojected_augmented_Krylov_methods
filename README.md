# Unprojected_augmented_Krylov_methods

A library which tests and compares unprojected augmented Krylov subspace methods with the corresponding projected augmented methods
and the standard (non augmented) method. 

This code is associated to the preprint 

[1] Burke, Liam, and Kirk M. Soodhalter. "Augmented unprojected Krylov subspace methods" arXiv preprint arXiv:2206.12315 (2022).

This library has the following file

* `fom.m`: A function which solves a linear system Ax = b using the standard restarted Full Orthogonalization Method (FOM).
* `r_fom.m`: A function which solves a linear system Ax = b using the (projected) recycled Full Orthogonalization Method (rFOM).
* `ur_fom.m`: A function which solves a linear system Ax = b using the (un projected) recycled Full Orthogonalization Method (urFOM).

* `gmres.m`: A function which solves a linear system Ax = b using the standard restarted Generalized Minimum Residual Method (GMRES).
* `r_gmres.m`: A function which solves a linear system Ax = b using the (projected) recycled Generalized Minimum Residual Method (rGMRES).
* `ur_gmres.m`: A function which solves a linear system Ax = b using the (un projected) recycled Generalized Minimum Residual Method (urGMRES).



## Tests

