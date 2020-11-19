#! /usr/bin/perl


# print stacked images
$usage = "usage: printstacked.pl ix1 ix2 iy1 iy2 dir prefix\
\tloop over ix = ix1 ... ix2; iy = iy1 ... iy\
\tstacked images should be in dir/ix_iy/prefix.fits\n";

$| = 1;

$filter = "imarith '%1 0.01 * 1.0 + fabs log 256 *' -";
$basedir = ".";


die $usage if ($#ARGV < 5);
$ix1 = shift @ARGV;
$ix2 = shift @ARGV;
$iy1 = shift @ARGV;
$iy2 = shift @ARGV;
$dir = shift @ARGV;
$prefix = shift @ARGV;

for ($ix = $ix1; $ix <= $ix2; $ix++) {
	for ($iy = $iy1; $iy <= $iy2; $iy++) {
		$X = ($ix - $offset) * $dX;
		$Y = ($iy - $offset) * $dY;
		$fitssrc = "$dir/$ix"."_$iy/$prefix.fits";
		$fitsdst = "/home/aaossl/nk/tmp/print$ix$iy.fits";
		$psfile = "/home/aaossl/nk/tmp/print$ix$iy.ps";
		&echosys("$filter < $fitssrc > $fitsdst");
		&echosys("print_image $fitsdst > $psfile");
		&echosys("lpr -s -h -Pnetps1 $psfile");
		&echosys("sleep 120");
		&echosys("rm $psfile $fitsdst");
	}
}


sub echosys {
        print @_, "\n";
        system(@_);
}
