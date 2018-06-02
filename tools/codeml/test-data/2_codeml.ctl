
      seqfile = /tmp/tmpr0j4bmyj/files/000/dataset_12.dat * sequence data file name
      outfile = run_codeml * main result file name
     treefile = /tmp/tmpr0j4bmyj/files/000/dataset_13.dat * tree structure file name
        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0  * 1: detailed output, 0: concise output
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise
      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 1  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0:no clock, 1:clock; 2:local clock
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
                   * 7:AAClasses
   aaRatefile = wag.dat * only used for aa seqs with model=empirical(_F)
                  * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
        model = 2
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional,2:Empirical,3:Empirical+F
                       * 6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)
      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0
        icode = 0  * 0:universal code; 1:mammalian mt; 2-11:see below
        Mgene = 0  * 0:rates, 1:separate;
    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2.0  * initial or fixed kappa
    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate
        omega = 0.2 * initial or fixed omega, for codons or codon-based AAs
    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0.0 * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * 1: different alphas for genes, 0 : one alpha
        ncatG = 3  * # of categories in dG of NSsites models
      fix_rho = 1  * 0: estimate rho; 1: fix it at rho
          rho = 0.0 * initial or fixed rho,   0:no correlation
        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
   Small_Diff = 5e-07
    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
  fix_blength = 0   * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 0   * 0: simultaneous; 1: one branch at a time
        