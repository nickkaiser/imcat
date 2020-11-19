#!/usr/bin/perl

# capture a set of frames and make a postscript file

$usage = "usage : $0 nx ny di dstps\n";

$nx 	= shift(@ARGV) || die $usage;
$ny 	= shift(@ARGV) || die $usage;
$di 	= shift(@ARGV) || die $usage;
$dstps 	= shift(@ARGV) || die $usage;

$basename = "/tmp/fits3Dviewer";
$i = 0;

$delta = 1.05;
$scale = 1.0 / ($delta * $nx);

$epsfcom = "epsfcompose -m 50 50 -a ";
for ($x = 0; $x < $nx; $x++) {
	for ($y = $ny - 1; $y >= 0; $y--) {
		$srcim = sprintf("%s.%.04d.ppm", $basename, $i);
		$psim = sprintf("%s.%.04d.ps", $basename, $i);
		echosys("pnmtops -noturn -rle < $srcim > $psim");
		$X = $x * $delta;
		$Y = $y * $delta;
		$epsfcom = $epsfcom."$psim $X $Y $scale 0 ";
		$i += $di;
	}
} 
echosys("$epsfcom > $dstps");

exit;

sub echosys {
	warn "$_[0]\n";
	system($_[0]) && die "system($_[0]) failed!\n";
}
