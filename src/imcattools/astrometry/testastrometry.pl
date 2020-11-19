#! /usr/bin/perl

# to make plots to check that the astrometry is working

# choose origin of coords
$ra0  = 0.0;
$dec0 = 0.0;

$ramin 	= -90.0;
$ramax 	=  90.0;
$Dra 	=  15.0;
$dra 	=   3.0;

$decmin = -75.0;
$decmax =  75.0;
$Ddec   =  15.0;
$ddec   =   3.0;

for ($ra = $ramin; $ra <= $ramax; $ra += $Dra) {
	for ($dec = $decmin; $dec <= $decmax; $dec += $ddec) {
		print `getxfromradec $ra $dec $ra0 $dec0`;
	}
}
for ($ra = $ramin; $ra <= $ramax; $ra += $dra) {
	for ($dec = $decmin; $dec <= $decmax; $dec += $Ddec) {
		print `getxfromradec $ra $dec $ra0 $dec0`;
	}
}