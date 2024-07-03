#!/usr/bin/env python3

import sys,numpy,math
from scipy import optimize
from scipy import special

Usage = """

extended_guinier.py N datafile outfile

where: N is the number of residues,
and data file has data points in the format q, I(q)

"""


def Rg_nu(b,N,nu):
    g = 1.1615
    return (g*(g+1.)/(2.*(g+2.*nu)*(g+2.*nu+1.)))**.5 * b * N**nu
    
def chisq_extguin(parms, q, log_Iq, Nres, b):
    nu = parms[0]
    log_I0 = parms[1]
    Rg = Rg_nu(b,Nres,nu)
    #rmax = 0.38*Nres

    log_Iq_fit = log_I0 - 1./3.*q**2*Rg**2 + 0.0479*(nu-0.212)*q**4 *Rg**4

    dI = log_Iq_fit - log_Iq
    chisq = numpy.dot(dI,dI)

    return chisq
    
b_Angstrom = 5.5 # fixed b for now

if len(sys.argv) != 4:
    sys.stdout.write( Usage )
    sys.exit(0)

Nres = int(sys.argv[1])
datafile = sys.argv[2]
outfile  = sys.argv[3]
data = numpy.loadtxt(datafile,comments=["#",])
q = data[:,0]
Iq = data[:,1]

#Rg_init = Rg_nu(b_Angstrom,Nres,0.5)
Rg_init = Rg_nu(b_Angstrom,Nres,0.6)
qcut = 2./Rg_init
q_lt_cut = q < qcut

qrange = q[q_lt_cut]
log_Iq = numpy.log(Iq[q_lt_cut])


nu_init = 0.5

log_I0_init = log_Iq[0]
init_parm = [nu_init,log_I0_init]

opt_parm = optimize.fmin(chisq_extguin, init_parm, (qrange,log_Iq,Nres,b_Angstrom))

nu_final,log_I0_final = opt_parm
Rg_final = Rg_nu(b_Angstrom,Nres,nu_final)
outp = open(outfile,"w")
outp.write("# =====================================\n")
outp.write("#                  nu = %8.3f\n"%(nu_final))
outp.write("# Extended Guinier Rg = %8.3f\n"%(Rg_final))
outp.write("# =====================================\n")
log_Iq_fit = log_I0_final - 1./3.*qrange**2*Rg_final**2 + 0.0479*(nu_final-0.212)*qrange**4 *Rg_final**4
for k,qk in enumerate(qrange):
    outp.write(f"{qk**2:12.6f} {log_Iq[k]:12.6f} {log_Iq_fit[k]:12.6f} \n")

outp.close()

