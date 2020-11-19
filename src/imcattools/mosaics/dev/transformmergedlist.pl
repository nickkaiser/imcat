#! /usr/bin/perl

# apply inverse transformations determined by mosaicfit to mergelists.tmp

require "chips.db";
require "fields.db";

@keys = keys %chipname;
$nchips = $#keys + 1;
@keys = keys %fieldname;
$nfields = $#keys + 1;


# read the transformation coefficients
open(COEFFTS, "mosaicfit.par") || die "Can't open mosaicfit.par";
($alpha) = split(' ', <COEFFTS>);
for ($fieldno = 0; $fieldno < $nfields; $fieldno++) {
	($Phi[$fieldno], $dX[$fieldno], $dY[$fieldno]) = split(' ', <COEFFTS>);
}
for ($chipno = 0; $chipno < $nchips; $chipno++) {
	($phi[$chipno], $dx[$chipno], $dy[$chipno]) = split(' ', <COEFFTS>);
}
close(COEFFTS);


# printf "%g\n", $alpha;
# for ($fieldno = 0; $fieldno < $nfields; $fieldno++) {
# 	printf "%g %g %g\n", $Phi[$fieldno], $dX[$fieldno], $dY[$fieldno];
# }
# for ($chipno = 0; $chipno < $nchips; $chipno++) {
# 	printf "%g %g %g\n", $phi[$chipno], $dx[$chipno], $dy[$chipno];
# }

while (<>) {
	($X[0], $Y[0], $c[0], $e[0], $X[1], $Y[1], $c[1], $e[1], $d) = split(' ');
	foreach $i (0,1) {
		$X[$i] -= $dX[$e[$i]];
		$Y[$i] -= $dY[$e[$i]];
		$Xe = $X[$i] - $Phi[$e[$i]] * $Y[$i];
		$Ye = $Y[$i] + $Phi[$e[$i]] * $X[$i];
		$rr = $Xe * $Xe + $Ye * $Ye;
		$xe = (1 - $alpha * $rr) * $Xe;
		$ye = (1 - $alpha * $rr) * $Ye;
		$xe -= $dx[$c[$i]];
		$ye -= $dy[$c[$i]];
		$X[$i] = $xe - $phi[$c[$i]] * $ye;
		$Y[$i] = $ye + $phi[$c[$i]] * $xe;
		printf "%13g %13g %2d %2d ", $X[$i], $Y[$i], $c[$i], $e[$i];
 	}
	printf "%13g\n", $d;
}


sub echosys {
	print @_, "\n";
	system(@_);
}
