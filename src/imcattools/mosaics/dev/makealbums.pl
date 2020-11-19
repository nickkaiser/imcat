#! /usr/bin/perl

# generate triply scrunched image of each exposure

$| = 1;

$N = 256;

require "chips.db";
require "fields.db";

$root = "fits";
# $root = "sky";
$scrunch = "scrunch | scrunch | scrunch";
# $scrunch = "scrunch";

# generate the offsets file
open(OFFSETS, ">tmp/album1.off") || die "Can't open tmp/album1.off\n";
printf OFFSETS "#	i	j\n";
foreach $c (keys %chipname) {
	printf OFFSETS "\t%d\t%d\n", 2 * $N * ($iy{$c} + 1), $N * ($ix{$c} + 2);
}
close(OFFSETS);
open(OFFSETS, ">tmp/album.off") || die "Can't open tmp/album.off\n";
printf OFFSETS "#\ti\tj\n";
printf OFFSETS "\t0\t0\n";
printf OFFSETS "\t$N\t0\n";
close(OFFSETS);

foreach $f (keys %fieldname) {
	$expname = $fieldname{$f};
	foreach $c (keys %chipname) {
		&echosys("cat < chip$chipname{$c}/$expname.$root | $scrunch > tmp/album.fits");
		if ($orient{$c} > 0) {
			&echosys("makesubimage 0 0 $N $N < tmp/album.fits > tmp/albuma$c.fits");
			&echosys("makesubimage $N 0 $N $N < tmp/album.fits > tmp/albumb$c.fits");
		} else {
			&echosys("makesubimage 0 0 $N $N < tmp/album.fits | spinflip -1 0 0 -1 > tmp/albuma$c.fits");
			&echosys("makesubimage $N 0 $N $N < tmp/album.fits | spinflip -1 0 0 -1 > tmp/albumb$c.fits");
		}
		if ($orient{$c} > 0) {
			&echosys("album -i 256 512 -f tmp/album.off tmp/albuma$c.fits tmp/albumb$c.fits > tmp/album$c.fits");
		} else {
			&echosys("album -i 256 512 -f tmp/album.off tmp/albumb$c.fits tmp/albuma$c.fits > tmp/album$c.fits");
		}
	}
	&echosys("album -i 1024 1024 -f tmp/album1.off -b tmp/album?.fits > albums/$expname.$root");
#	&echosys("print_image albums/$expname.$root -f -30 200 > albums/$expname.$root.ps");
#	&echosys("lpr -h -s -Pnetps1 albums/$expname.$root.ps");
}
 	
&echosys("rm tmp/album*");	


sub echosys {
        print @_, "\n";
        system(@_);
}
