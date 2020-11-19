#! /usr/bin/perl

$ncom = 4;

echosys("ic -c 32 32 'grand' > tmp0.fits");

$command = "ic -h BSCALE 2.3 ";

while ($ncom > 0) {
	$ncom--;
	$command = sprintf("%s -h FOO_$ncom BAR_$ncom", $command);
}
$command = sprintf("%s '%%1' tmp0.fits > tmp.fits", $command);

echosys($command);
echosys("rm tmp0.fits");



sub echosys() {
	print @_[0], "\n";
	system(@_[0]) && die "failed to execute @_[0]\n";
}
