#!/usr/bin/perl

$dt = 3.14159 / 180.0;
$sigma = 5.e-7;
$dtmax = 0.5 * $dt;

$np = 256;

# cadence in days
$nsteps = 4;

# number of observations
$nt = 4.0;

$Dt = $nsteps * $dt;

$srccat = "tmp_fitorbit_src.cat";
$initsolcat = "tmp_fitorbit_ini.cat";
$obscat = "tmp_fitorbit_obs.cat";
$outcat = "tmp_fitorbit_out.cat";

chop($fnummin = `lc fnum < $srccat | catstats -v min`);

if (0) {

# chop($nt = `lc -c < $srccat`);

# generate the initial solution catalog
echosys("lc -b -i '%fnum $fnummin ==' 'ra0 = %ra' 'va0 = %va' 're0 = %re' 've0 = %ve' 't0 = 0.0' 'nt = $nt' < $srccat > $initsolcat");

# evolve the earth and the asteroid
foreach $tag ("a", "e") {
	echosys("lc -b -i '%fnum $fnummin ==' 'r = %r$tag' 'v = %v$tag' 'risk = 1.0' < $srccat ".
		"| ./tcl_evolveN $dt $nsteps $nt | lc -b 'r$tag = %r' 'v$tag = %v' +cI > tmp.$tag.cat");
}

# generate the observations catalog
echosys("pastecats tmp.[ae].cat | lc -b 't = %I $Dt *' re 'rho = 0 0 0 3 vector' ".
	"'nobs = %ra %re vsub enter enter dot -0.5 pow vscale' ".
	"'sigma = $sigma' > $obscat"); 

echosys("./fitorbit $initsolcat $dtmax $np < $obscat | lc -b +all 'risk = 1' > $outcat");

}

echosys("./tcl_evolveN $dt 720 1 < $outcat | ./orbs2obs $np 1 $dt > tmp.orbs2obs.cat");

sub echosys {
        warn "$_[0]\n";
        system($_[0]) && die "system($_[0]) failed\n";
}
