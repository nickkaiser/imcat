#! /usr/bin/perl

# apply transformations defined by mosaicfit.par to a superlist

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


while (<>) {
	unless (/^#/) {
		($xce, $yce, $c, $e) = split(' ');
		$xe = $xce + $phi[$c] * $yce + $dx[$c];
		$ye = $yce - $phi[$c] * $xce + $dy[$c];
		$rr = $xe * $xe + $ye * $ye;
		$Xe = (1 + $alpha * $rr) * $xe;
		$Ye = (1 + $alpha * $rr) * $ye;
		$X = $Xe + $Phi[$e] * $Ye + $dX[$e];
		$Y = $Ye - $Phi[$e] * $Xe + $dY[$e];
		printf "%13g %13g %d %d\n", $X, $Y, $c, $e;
	}
}
