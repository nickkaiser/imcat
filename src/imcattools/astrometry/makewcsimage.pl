#!/usr/bin/perl

$PI = 3.14158265358979323846;
$rad = 180.0 / $PI;
$type = "-STG";
$ra0 = 0.0;
$dec0 = 0.0;
$N = 240;
$delt = 1.0;

$N2 = $N / 2;

echosys("ic -c $N $N ".
	"-h CTYPE1 RA--$type ".
	"-h CDELT1 $delt ".
	"-h CRVAL1 $ra0 ".
	"-h CRPIX1 $N2 ".
	"-h CTYPE2 DEC-$type ".
	"-h CDELT2 $delt ".
	"-h CRVAL2 $dec0 ".
	"-h CRPIX2 $N2 ".
	"'grand' > tmp.fits");
 

sub echosys {
        print @_, "\n";
        system(@_);
}
