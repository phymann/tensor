# tensor

Examples of using methods like DMRG, MPS, Tensor network, etc, that involve manipulating tensors.



## Gallery

1. [TRG.jl](https://github.com/phymann/tensor/blob/master/TRG.jl)
	
	​		Free energy per site for classical 2D Ising model is obtained using tensor renormalization group algorithm proposed in [cond-mat/0611687](http://arxiv.org/abs/cond-mat/0611687). This code is essentially a rewritten of [C++ version](http://itensor.org/docs.cgi?vers=cppv3&page=book/trg) into a Julia version.

2. [TwoCompBHM.ipynb](https://github.com/phymann/tensor/blob/master/TwoCompBHM.ipynb)

   ​		This IJulia notebook works on the two-component Bose-Hubbard model in the strong coupling limit at total filling $2N_{\text{site}}$, whose effective Hamiltonian is a modified ferro-Heisenberg spin model (the derivation is given [here](https://github.com/phymann/Miscellaneous/blob/master/effective_hamiltonian.nb)). Details about this problem can be found in these two papers:

   1. [Phase diagram of two-component bosons on an optical lattice](https://iopscience.iop.org/article/10.1088/1367-2630/5/1/113)
   2. [Adiabatic cooling of bosons in lattices to magnetically ordered quantum states](https://link.aps.org/doi/10.1103/PhysRevA.92.041602)

