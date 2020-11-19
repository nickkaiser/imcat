#!/usr/bin/perl

# exercise testkepler
$usage = "usage : $0 norbits\n";

$norbits = shift(@ARGV) || die $usage;

echosys("makerandcat $norbits -dim 6 | lc 'x = %x -0.5 vshift 0.05 vscale' ".
	"| lc 'r = %x[0] 1 + %x[1] %x[2] 3 vector' 'v = %x[3] %x[4] 1 + %x[5] 3 vector' ".
#	"| ./testkepler | lc r_in r_out v_in v_out");
	"| ./testkepler");

sub echosys {
	warn "$_[0]\n";
	system($_[0]) && die "system($_[0]) failed\n";
}
