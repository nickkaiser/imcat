#! /usr/bin/perl

# generate a stack of images for a given region of X,Y sky coords
$usage = "usage: makestack.pl X Y dX dY destdir\
\twarp images to create a stack for sub-image of size dX, dY at X, Y\n";

$| = 1;

$basedir = ".";
$destdir = "stacked";


die $usage if ($#ARGV < 4);
$X[0] = shift @ARGV;
$Y[0] = shift @ARGV;
$dX = shift @ARGV;
$dY = shift @ARGV;
die "non-positive sub-image size\n" if ($dX <= 0 || $dY <= 0);
$X[1] = $X[0] + $dX;
$Y[1] = $Y[0] + $dX;
$destdir = shift @ARGV;

require "chips.db";
require "fields.db";
require "nominal.db";

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

foreach $f (keys %fieldname) {
	$target = "$destdir/$fieldname{$f}.fits";
	&echosys("make_image 1 $dX $dY | remagic -f -1 1 > $target");
	foreach $c (keys %chipname) {
		$use[$c] = 0;
	}
	foreach $c (keys %chipname) {
		foreach $h (0,1) {
			foreach $v (0,1) {
				$x = $X[$h] - $dX[$f];
				$y = $Y[$v] - $dY[$f];
				$xx = $x - $Phi[$f] * $y;
				$yy = $y + $Phi[$f] * $x;
				$rr = $xx * $xx + $yy * $yy;
				$x = (1 - $alpha * $rr) * $xx;
				$y = (1 - $alpha * $rr) * $yy;
				$x -= $dx[$c];
				$y -= $dy[$c];
				$xx = $x - $phi[$c] * $y;
				$yy = $y + $phi[$c] * $x;
				$x = $xx - $ix{$c} * $Xsize - $xmargin;
				$y = $yy - $iy{$c} * $Ysize - $ymargin;
				if (($x >= 0.0 && $x < $xsize) && ($y >= 0.0 && $y < $ysize)) {
					$use[$c]++;
				}
			}
		}
	}
	foreach $c (keys %chipname) {
		if ($use[$c] > 0) {
			$src = "$basedir/chip$chipname{$c}/$fieldname{$f}.fits";
			$sky = "$basedir/chip$chipname{$c}/$fieldname{$f}.sky";
			$DX = $ix{$c} * $Xsize + $xmargin;
			$DY = $iy{$c} * $Ysize + $ymargin;
			&echosys("unscrunch < $sky | unscrunch | ic '%1 %2 -' $src - > $destdir/temp.fits");
			&echosys("mosaicmap  $alpha $Phi[$f] $dX[$f] $dY[$f] $phi[$c] $dx[$c] $dy[$c] ".
				"$DX $DY $X[0] $Y[0] $orient{$c} $xsize $ysize $target < $destdir/temp.fits");
		}
	}
}
&echosys("rm $destdir/temp.fits");
&echosys("pastiche -i $dX $dY -c $destdir/0??.fits > $destdir/count.fits");
&echosys("pastiche -i $dX $dY $destdir/0??.fits > $destdir/median.fits");
&echosys("pastiche -i $dX $dY -a 3 $destdir/0??.fits > $destdir/avsgclp.fits");
&echosys("rm $destdir/0??.fits");



sub echosys {
        print @_, "\n";
        system(@_);
}
