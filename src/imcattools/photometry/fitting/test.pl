#! /usr/bin/perl

# number of objects
$nobj = 3;

# image size 
$N1 = 18;
$N2 = 24;

# maximum central value
$f0max = 3000.0;

# axial ratios
$rmin = 3;
$rmax = 5;

srand(100);

open(CAT, "|lc -C -n f0 -N '1 2 x' -n a -n b -n phi > test.cat");
echosys("ic -c $N1 $N2 '0' > test.fits");
while ($nobj--) {
	$f0 = $f0max * rand();
	$x = $N1 * rand();
	$y = $N2 * rand();
	$a = $rmin + ($rmax - $rmin) * rand();
	$b = $rmin + ($rmax - $rmin) * rand();
	$phi = 3.14159 * rand();
	print CAT "$f0 $x $y $a $b $phi\n";
	$a2 = $a * $a;
	$b2 = $b * $b;
	$c = cos($phi);
	$s = sin($phi);
	echosys("ic -c $N1 $N2 'xp $x - $c * yp $y - $s * + enter * $a2 / yp $y - $c * xp $x - $s * - enter * $b2 / + -0.5 * exp $f0 *' > tmp.fits");
	echosys("doto test.fits ic \"'%1 %2 +'\" tmp.fits -");
}
close(CAT);
echosys("rm tmp.fits");

sub echosys {
        print @_, "\n";
        system(@_) && die "system call \"@_\" failed\n";
}
