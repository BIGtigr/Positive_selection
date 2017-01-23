#!/usr/bin/perl

use strict;
use warnings;

my($Alignment) = $ARGV[0];
my($codeml_output) = $ARGV[1];
my($model) = $ARGV[2];
my($tree) = $ARGV[3];

## HASH model, NSsites, ncatG, fix_omega, omega
my(%models) = (
	'Model01ratio' 		=> [ 0, 0, 10, 0, 0.4],   ################ if wanting to run EVOLVER on the output remember to put "verbose = 2"
	'Model1Neutral' 	=> [ 0, 1, 10, 0, 0.4],
	'Model2Selection' 	=> [ 0, 2, 10, 0, 0.4],
	'Model3Discrtk2' 	=> [ 0, 3, 2, 0, 0.4],
	'Model3Discrtk3' 	=> [ 0, 3, 3, 0, 0.4],
	'Model7beta' 		=> [ 0, 7, 10, 0, 0.4],
	'Model8beta' 		=> [ 0, 8, 10, 0, 0.4],
	'ModelAbranchSite' 	=> [ 2, 2, 3, 0, 0.4],
	'NullModelAbranchSite' 	=> [ 2, 2, 3, 1, 1],
	'ModelBbranchSite' 	=> [ 2, 3, 3, 0, 0.4],
	'CladeModelC'		=> [ 3, 2, 3, 0, 0.4],
	'NullCladeModelC'		=> [ 3, 2, 3, 1, 1],
	'Model2RatioBranch'	=> [ 2, 0, 10, 0, 0.4],
	'ModelFreeRatio' 	=> [ 1, 0, 10, 0, 0.4],
		);

#open(OUT,">codeml.ctl") || die;

if(exists $models{$model}){
	open(OUT,">codeml.ctl") || die;
	print OUT "          seqfile  = $Alignment
          treefile = $tree
          outfile  = $codeml_output

             noisy = 0  * 0,1,2,3,9: how much rubbish on the screen
           verbose = 0  * 0: concise; 1: detailed, 2: too much
           runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                        * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

           seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
         CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
            aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a

        aaRatefile = wag.dat * only used for aa seqs with model=empirical(_F)
                             * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
             model = $models{$model}[0] * models for codons:
                              * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                              * models for AAs or codon-translated AAs:
                              * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                              * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)
           NSsites = $models{$model}[1] * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                              * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                              * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                              * 13:3normal>0

             icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
         fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
             kappa = 3  * initial or fixed kappa
         fix_omega = $models{$model}[3]  * 1: omega or omega_1 fixed, 0: estimate
             omega = $models{$model}[4] * initial or fixed omega, for codons or codon-based AAs

         fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
             alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
            Malpha = 0  * different alphas for genes

             ncatG = $models{$model}[2] * # of categories in dG of NSsites models
             clock = 0  * 0:no clock, 1:global clock; 2:local clock; 3:TipDate
             getSE = 0  * 0: don't want them, 1: want S.E.s of estimates

      RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
        Small_Diff = .5e-6\n";

	close OUT;

	system("codeml");

}elsif($model eq "ModelPairwise"){
	open(OUT,">codeml.ctl") || die;
	print OUT "          seqfile  = $Alignment
          treefile = 
          outfile  = $codeml_output

             noisy = 0  * 0,1,2,3,9: how much rubbish on the screen
           verbose = 0  * 0: concise; 1: detailed, 2: too much
           runmode = -2  * 0: user tree;  1: semi-automatic;  2: automatic
                        * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

           seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
         CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
            aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a

        aaRatefile = wag.dat * only used for aa seqs with model=empirical(_F)
                             * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
             model = 0  * models for codons:
                              * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                              * models for AAs or codon-translated AAs:
                              * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                              * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)
           NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                              * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                              * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                              * 13:3normal>0

             icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
         fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
             kappa = 2  * initial or fixed kappa
         fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate
             omega = 0.4 * initial or fixed omega, for codons or codon-based AAs

         fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
             alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
            Malpha = 0  * different alphas for genes

             ncatG = 10 * # of categories in dG of NSsites models
             clock = 0  * 0:no clock, 1:global clock; 2:local clock; 3:TipDate
             getSE = 0  * 0: don't want them, 1: want S.E.s of estimates

      RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
        Small_Diff = .5e-6\n";

	close OUT;

	system("codeml");

}elsif($model eq "YN00"){
	open(OUT,">yn00.ctl") || die;
	print OUT "           seqfile = $Alignment
           outfile = $codeml_output
           verbose = 1
             icode = 0
         weighting = 0
        commonf3x4 = 0";

	close OUT;

	system("yn00");
}

#close OUT;

#system("codeml");

print "finished\n";

exit;
