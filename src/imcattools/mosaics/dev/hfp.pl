#! /usr/bin/perl


# run hfindpeaks on stacked images
$usage = "usage: hfp.pl ix1 ix2 iy1 iy2 dir prefix\
\tloop over ix = ix1 ... ix2; iy = iy1 ... iy\
\tstacked images should be in dir/ix_iy/prefix.fits\n";

$| = 1;			# force flushing



die $usage if ($#ARGV < 5);
$ix1 = shift @ARGV;
$ix2 = shift @ARGV;
$iy1 = shift @ARGV;
$iy2 = shift @ARGV;
$dir = shift @ARGV;
$prefix = shift @ARGV;

for ($ix = $ix1; $ix <= $ix2; $ix++) {
	for ($iy = $iy1; $iy <= $iy2; $iy++) {
		$fitsfile = "$dir/$ix"."_$iy/$prefix.fits";
		$catfile = "$dir/$ix"."_$iy/$prefix.cat";
		&echosys("hfindpeaks $fitsfile > $catfile");
	}
}


sub echosys {
        print @_, "\n";
        system(@_);
}
