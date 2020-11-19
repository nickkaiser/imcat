#! /usr/bin/perl


$offset = 5;


# generate a grid of stacked images
$usage = "usage: makestack.pl ix1 ix2 iy1 iy2 dX dY destdir\
\tgenerate a grid of stacked images of size dX, dY.\
\tloop over ix = ix1 ... ix2; iy = iy1 ... iy\
\torigin is at (ix - $offset) * dX, (iy - $offset) * dY\
\tstacked images placed in destdir/ix_iy/\n";

$| = 1;

$basedir = ".";


die $usage if ($#ARGV < 6);
$ix1 = shift @ARGV;
$ix2 = shift @ARGV;
$iy1 = shift @ARGV;
$iy2 = shift @ARGV;
$dX = shift @ARGV;
$dY = shift @ARGV;
die "non-positive sub-image size\n" if ($dX <= 0 || $dY <= 0);
$destdir = shift @ARGV;

for ($ix = $ix1; $ix <= $ix2; $ix++) {
	for ($iy = $iy1; $iy <= $iy2; $iy++) {
		$X = ($ix - $offset) * $dX;
		$Y = ($iy - $offset) * $dY;
		$subdir = "$destdir/$ix"."_$iy";
		&echosys("mkdir $subdir");
		&echosys("makestack.pl $X $Y $dX $dY $subdir");
	}
}


sub echosys {
        print @_, "\n";
       system(@_);
}
