## 2D pore-network model

A pore-network model is a geometrical approximnation of complex porous media where large pores are defined as pore-bodies to store fluids and narrow channels between them are defined as pore-throats.

A two-dimensional pore-network model generator with arbitrary pore-body and pore-throat size distribution is written in MATLAB. The code calculates single-phase flow through pore-throats and solves for pressure in pore-bodies. The ouput is a representation of the generated pore-network with predicted pressure field and its average distribution.


More details on pore-network modeling can be found in these sources:

[1] Fatt, I. (1956). [The network model of porous media](https://www.onepetro.org/general/SPE-574-G). Petroleum Transactions, AIME, 207, 144-181.

[2] Blunt, M. J. (2001). [Flow in porous media—pore-network models and multiphase flow](https://doi.org/10.1016/S1359-0294(01)00084-X). Current opinion in colloid & interface science, 6(3), 197-207.

[3] Valvatne, P. H., & Blunt, M. J. (2004). [Predictive pore‐scale modeling of two‐phase flow in mixed wet media](https://doi.org/10.1029/2003WR002627). Water resources research, 40(7).

[4] Raeini, A. Q., Bijeljic, B., & Blunt, M. J. (2017). [Generalized network modeling: Network extraction as a coarse-scale discretization of the void space of porous media](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.96.013312). Physical Review E, 96(1), 013312.
