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
from ActiveSpaceMaker import MakeActiveSpace
from pyscf import mrpt



def readgeom(xyzfile):
 input=open(xyzfile,'r')
 lines=input.readlines()
 input.close()
 N_lines=len(lines)
 AtomsPyscf = []
 for n in range (2,N_lines):
  AtomsPyscf.append(lines[n])
 AtomsPyscf = "\n".join(AtomsPyscf)
 return AtomsPyscf

mol = gto.Mole()
xyz_input= "bluf_oh_10_nh_10.xyz"
mol.atom=readgeom(xyz_input)
mol.verbose =4
mol.basis = '6-31++G**'
mol.spin = 0 #na-nb
mol.build()

mf=scf.ROHF(mol)
mf.chkfile = 'hf_oh_1.0_nh_1.0.chk'
mf.kernel()

N_Core,N_Act, N_Virt,nelec,coeff=MakeActiveSpace(mol,mf)
print " # of core orbs    ", N_Core
print " # of active orbs  ", N_Act
print " # of virtual orbs ", N_Virt
nalpha=(nelec+mol.spin)/2
nbeta=(nelec-mol.spin)/2


mycas = mcscf.CASSCF(mf, N_Act, [nalpha,nbeta])
AS=range(N_Core,N_Core+N_Act)
mycas =mycas.state_average_([0.167,0.1666,0.1666,0.1666,0.1666,0.1666])
mycas.chkfile ='cas_oh_1.0_nh_1.0.chk'
mycas.fcisolver.nroots = 6
mycas.fix_spin_(ss=0)
activeMO = mcscf.sort_mo(mycas,coeff,AS,base=0)
mycas.verbose = 5
mycas.max_cycle_macro = 150
mycas.kernel(activeMO)

mc = mcscf.CASCI(mf, N_Act, [nalpha,nbeta])
mc.__dict__.update(scf.chkfile.load('cas_oh_1.0_nh_1.0.chk', 'mcscf'))
mc.fcisolver.nroots = 6
mc.fix_spin_(ss=0)
mc.verbose = 5
mc.kernel()


def analyze_charges(Charges,n_excited_states):
 # here for charge analysis I did not includ H atoms, since the main redistribution of charges happens on pi-system C, N, O atoms
   TrpPiAtoms = [27,28,30,32,33,34,36,38,40]
   # FmnPiAtoms = [17, 18, 20, 21, 23, 24, 25, 26, 27, 29, 31, 32, 33, 34]
   # note O 19 and O 23 not part of ring, but their lone pairs are clearly delocalized into ring.
   FmnPiAtoms = [42,43,44,45,47, 48,49,50,51,52, 54,59,64,66,67, 68]
   TyrfPiAtoms = [4,5,7,9,12,14]
   GlnAtoms = [16,19,20,21]

   chg_Trp= numpy.zeros(5)
   chg_Fmn=numpy.zeros(5)
   chg_Tyr=numpy.zeros(5)
   chg_Gln=numpy.zeros(5)

   n_Trp=len(TrpPiAtoms)
   n_Tyr=len(TyrPiAtoms)
   n_Fmn=len(FmnPiAtoms)
   n_Gln=len(GlnPiAtoms)

   for state in range(1,n_excited_states+1):
    dchg_Trp=0.0
    dchg_Tyr=0.0
    dchg_Fmn=0.0
    dchg_Gln=0.0

    for k in range(n_Trp):
     ia=TrpPiAtoms[k]-1
     dchg_Trp+=(Charges[state][ia]-Charges[0][ia])
    chg_Trp.append(dchg_Trp)

    for k in range(n_Tyr):
     ia=TyrPiAtoms[k]-1
     dchg_Tyr+=(Charges[state][ia]-Charges[0][ia])
    chg_Tyr.append(dchg_Tyr)

    for k in range(n_Fmn):
     ia=FmnPiAtoms[k]-1
     dchg_Fmn+=(Charges[state][ia]-Charges[0][ia])
    chg_Fmn.append(dchg_Fmn)

    for k in range(n_Gln):
     ia=GlnPiAtoms[k]-1
     dchg_Gln+=(Charges[state][ia]-Charges[0][ia])
    chg_Gln.append(dchg_Gln)

   for i in range(n_excited_states):
      print ('relative charges for Fmn, Tyr, Trp and Gln residues for excited state # %d (w.r.t. the Ground State) are %5.5f    %5.5f   %5.5f   %5.5f  '% (i+1,chg_Fmn[i],chg_Tyr[i],chg_Trp[i],chg_Gln[i]))

def read_charges(casci_object,n_excited_states):
  chgs=[]
  for i in range(1+n_excited_states):
    ncore = mc.ncore
    ncas = mc.ncas
    nocc = ncore + ncas
    nelecas = mc.nelecas
    casdm1= mc.fcisolver.make_rdm1(mc.ci[i], ncas, nelecas)
    mo_coeff = mc.mo_coeff
    mocore = mo_coeff[:,:ncore]
    mocas = mo_coeff[:,ncore:nocc]
    dm1 =(numpy.dot(mocore, mocore.T) * 2
             + reduce(numpy.dot, (mocas, casdm1, mocas.T)))
    ovlp_ao = mc._scf.get_ovlp()
    pop,chg= mc._scf.mulliken_meta(mc.mol, dm1, s=ovlp_ao)
    chgs.append(chg)
  return chgs

charges=read_charges(mc,5)
analyze_charges(charges,5)

# compute NEVPT2 corrections
ci_nevpt_e1 = mrpt.NEVPT(mc, root=0).kernel()
ci_nevpt_e2 = mrpt.NEVPT(mc, root=1).kernel()
ci_nevpt_e3 = mrpt.NEVPT(mc, root=2).kernel()
ci_nevpt_e4 = mrpt.NEVPT(mc, root=3).kernel()
ci_nevpt_e5 = mrpt.NEVPT(mc, root=4).kernel()
ci_nevpt_e6 = mrpt.NEVPT(mc, root=5).kernel()



