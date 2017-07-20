#!/usr/bin/env python

import os, sys

# Command line : python S01_comp_codeml.py <tree_file> <concat_file> <model> <NSsites> <fix_omega> <omega>

tree_file = sys.argv[1]
concat_file = sys.argv[2]
model = sys.argv[3]
NSsites = sys.argv[4]
fix_omega = sys.argv[5]
omega = sys.argv[6]

print "tree : %s" %tree_file
print "sequences : %s" %concat_file
print "model : %s" %model
print "NSsites : %s" %NSsites
print "fix_omega : %s" %fix_omega
print "omega : %s" %omega

codeml=open("codeml.ctl", "w")
codeml.write("      seqfile = %s * sequence data file name \n" %concat_file )
codeml.write("      outfile = run_codeml * main result file name \n")
codeml.write("     treefile = %s * tree structure file name\n" %tree_file)
codeml.write("        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen \n")
codeml.write("      verbose = 0  * 1: detailed output, 0: concise output \n")
codeml.write("      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic \n")
codeml.write("                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise \n")
codeml.write("      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs \n")
codeml.write("    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table \n")
codeml.write("*       ndata = 10 \n")
codeml.write("        clock = 0   * 0:no clock, 1:clock; 2:local clock \n")
codeml.write("       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a \n")
codeml.write("                   * 7:AAClasses \n")
codeml.write("   aaRatefile = wag.dat * only used for aa seqs with model=empirical(_F) \n")
codeml.write("                  * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own \n")
codeml.write("        model = %s \n" %model)
codeml.write("                   * models for codons: \n")
codeml.write("                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches \n")
codeml.write("                   * models for AAs or codon-translated AAs: \n")
codeml.write("                       * 0:poisson, 1:proportional,2:Empirical,3:Empirical+F \n")
codeml.write("                       * 6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189) \n")
codeml.write("      NSsites = %s  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs; \n" %NSsites)
codeml.write("                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma; \n")
codeml.write("                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1; \n")
codeml.write("                   * 13:3normal>0 \n")
codeml.write("        icode = 0  * 0:universal code; 1:mammalian mt; 2-11:see below \n")
codeml.write("        Mgene = 0  * 0:rates, 1:separate;  \n")
codeml.write("    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated \n") 
codeml.write("        kappa = 2  * initial or fixed kappa\n")
codeml.write("    fix_omega = %s  * 1: omega or omega_1 fixed, 0: estimate  \n" %fix_omega)
codeml.write("        omega = %s * initial or fixed omega, for codons or codon-based AAs \n" %omega)
codeml.write("    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha \n")
codeml.write("        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate) \n")
codeml.write("       Malpha = 0  * different alphas for genes \n")
codeml.write("        ncatG = 3  * # of categories in dG of NSsites models \n")
codeml.write("      fix_rho = 1  * 0: estimate rho; 1: fix it at rho \n")
codeml.write("          rho = 0. * initial or fixed rho,   0:no correlation \n")
codeml.write("        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates \n")
codeml.write(" RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2) \n")
codeml.write("   Small_Diff = .5e-6 \n")
codeml.write("*   cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)? \n")
codeml.write("* fix_blength = 0   * 0: ignore, -1: random, 1: initial, 2: fixed \n")
codeml.write("       method = 0   * 0: simultaneous; 1: one branch at a time ")
codeml.close()

os.system("codeml codeml.ctl")
