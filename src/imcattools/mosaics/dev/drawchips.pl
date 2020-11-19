#! /usr/bin/perl

require 'read_config_db.pl';

$margin = 30;
$xsize = 2048;
$ysize = 4096;
$Xsize = $xsize + 2 * $margin;
$Ysize = $ysize + 2 * $margin;

&read_config_db;



foreach $thechipno (keys %ix) {
	$cos = $cos{$thechipno};
	$x = $Xsize * $ix{$thechipno} + $cos * $margin;
	$y = $Ysize * $iy{$thechipno} + $cos * $margin;
	printf "rel %g %g\ndraw %g %g\ndraw %g %g\ndraw %g %g\ndraw %g %g\n", 
		$x, $y, 
		$x + $cos * $xsize, $y,
		$x + $cos * $xsize, $y + $cos * $ysize,
		$x, $y + $cos * $ysize,
		$x, $y;
}


