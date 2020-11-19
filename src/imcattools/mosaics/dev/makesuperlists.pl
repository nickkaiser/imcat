#! /usr/bin/perl

# create superlists

$catdir = ".";
# $catdir = "mock";
$superlistdir = "superlists";

require "nominal.db";
require "chips.db";
require "fields.db";

foreach $fieldnum (keys %fieldname) {
	print "# field $fieldnum...\n";
	$thesuperlist = "$superlistdir/field$fieldnum.sl";
	open(SL, ">$thesuperlist") || die "Can't open $thesuperlist\n";
	foreach $chipnum (keys %chipname) {
		$dx = $ix{$chipnum} * $Xsize + $xmargin;
		$dy = $iy{$chipnum} * $Ysize + $ymargin;
		$thecat = "$catdir/chip$chipname{$chipnum}/$fieldname{$fieldnum}.stars";
		open(CAT, "listcat < $thecat |") || die "Can't open $thecat\n";
		while (<CAT>) {
			unless (/^#/) {
				($y,$x) = split(' ');
				if ($orient{$chipnum} < 1) {
					$x = $xsize - 1 - $x;
					$y = $ysize - 1 - $y;
				}
				$x += $dx;
				$y += $dy;
				printf SL "%13g %13g %2d %2d\n", 
						$x, $y, $chipnum, $fieldnum;
			}
		}
		close(CAT);
	}
	close(SL);
}

exit;





sub echosys {
	print @_, "\n";
	system(@_);
}
