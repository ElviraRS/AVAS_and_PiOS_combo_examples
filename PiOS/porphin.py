import numpy
from pyscf import gto
from pyscf.gto import mole
from pyscf.gto import moleintor
from pyscf import scf
from pyscf import ao2mo
from pyscf import mcscf
from functools import reduce
from pyscf import fci
from pyscf.tools import fcidump
from pyscf import lib
from pyscf import mrpt
from PiOS import MakePiOS


mol = gto.Mole()
mol.atom = '''
C -0.682956 -4.256282 0.000000
C -1.083580 -2.853519 0.000000
N 0.000000 -2.024853 0.000000
C 1.083580 -2.853519 0.000000
C 0.682956 -4.256282 0.000000
C -2.425952 -2.443542 0.000000
C -2.900522 -1.126914 0.000000
C -4.268346 -0.689607 0.000000
C -4.268346 0.689607 0.000000
C -2.900522 1.126914 0.000000
N -2.123164 0.000000 0.000000
C -2.425952 2.443542 0.000000
C -1.083580 2.853519 0.000000
C -0.682956 4.256282 0.000000
C 0.682956 4.256282 0.000000
C 1.083580 2.853519 0.000000
N 0.000000 2.024853 0.000000
C 2.425952 2.443542 0.000000
C 2.900522 1.126914 0.000000
C 4.268346 0.689607 0.000000
C 4.268346 -0.689607 0.000000
C 2.900522 -1.126914 0.000000
N 2.123164 0.000000 0.000000
C 2.425952 -2.443542 0.000000
H 5.133687 -1.357782 0.000000
H 5.133687 1.357782 0.000000
H 1.362163 5.114185 0.000000
H -1.362163 5.114185 0.000000
H -5.133687 1.357782 0.000000
H -5.133687 -1.357782 0.000000
H -1.362163 -5.114185 0.000000
H 1.362163 -5.114185 0.000000
H 3.190898 -3.230970 0.000000
H 3.190898 3.230970 0.000000
H -3.190898 3.230970 0.000000
H -3.190898 -3.230970 0.000000
H -1.096774 0.000000 0.000000
H 1.096774 0.000000 0.000000
'''
mol.basis = 'cc-pvdz'
mol.spin = 0
mol.build()

mf = scf.RHF(mol)
mf.chkfile = 'hf_porphin.chk'
mf.kernel()

PiAtoms=range(1,25)
N_Core,N_Act, N_Virt,nelec,coeff=N_Core,N_Act, N_Virt,nelec,coeff=MakePiOS(mol,mf,PiAtoms,2,2)
print " # of core orbs    ", N_Core
print " # of active orbs  ", N_Act
print " # of virtual orbs ", N_Virt
nalpha=(nelec+mol.spin)/2
nbeta=(nelec-mol.spin)/2


mycas = mcscf.CASSCF(mf, N_Act, [nalpha,nbeta])
AS=range(N_Core,N_Core+N_Act)
mycas =mycas.state_average_([0.2,0.2,0.2,0.2,0.2])
mycas.chkfile ='cas_porphin.chk'
mycas.fcisolver.nroots = 5
activeMO = mcscf.sort_mo(mycas,coeff,AS,base=0)
mycas.fix_spin_(ss=0)
mycas.verbose = 5
mycas.max_cycle_macro = 30
mycas.kernel(activeMO)

mycas = mcscf.CASCI(mf, N_Act, [nalpha,nbeta])
mycas.__dict__.update(scf.chkfile.load('cas_porphin.chk', 'mcscf'))
mycas.fix_spin_(ss=0)
mycas.fcisolver.nroots = 5
mycas.verbose = 5
mycas.kernel()

ci_nevpt_e1 = mrpt.NEVPT(mycas, root=0).kernel()
ci_nevpt_e2 = mrpt.NEVPT(mycas, root=1).kernel()
ci_nevpt_e3 = mrpt.NEVPT(mycas, root=2).kernel()
ci_nevpt_e4 = mrpt.NEVPT(mycas, root=3).kernel()
ci_nevpt_e5 = mrpt.NEVPT(mycas, root=4).kernel()
