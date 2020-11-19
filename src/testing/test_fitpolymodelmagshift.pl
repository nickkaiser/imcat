#!/usr/bin/perl

$nobj = 10000;
$mag0 = 15;
# $magrange = 5;
$magrange = 0;
$imsize = 10000.0;
# $imsize = 1.0;
$xscale = 2 * $imsize;
$fitorder = 4;
$nrand = 100;
$deltamagscale = 0.1;
$halfimsize = 0.5 * $imsize;
$xydither = 0.2 * $imsize;
$mergetol = 1.e-3 * $imsize;

$flaterror = "| lc -b +all 'mag = %mag %xdet %xdet dot $halfimsize enter * / $deltamagscale * +' ";

# generate master catalogue
echosys("makerandcat $nobj | lc -b 'x = %x -0.5 vshift $xscale vscale' 'mag = rand $magrange * $mag0 +' > foo.master.cat");
# make the subsets
$e = 0;
echosys("lc -b -i '%x[0] fabs $halfimsize < %x[1] fabs $halfimsize < and' +all 'xdet = %x' < foo.master.cat $flaterror 'e = $e' > foo.$e.cat");
$e++;
foreach $scalefac (1) {
	foreach $iy (-1, 1) {
		$dy = $scalefac * $xydither * $iy;
		foreach $ix (-1, 1) {
			$dx = $scalefac * $xydither * $ix;
			echosys("lc -b -i '%x[0] $dx + fabs $halfimsize < %x[1] $dy + fabs $halfimsize < and' < foo.master.cat ".
				"| lc -b +all 'xdet = %x $dx $dy 2 vector vadd' $flaterror 'e = $e' > foo.$e.cat");
			$e++;
		}
	}
}
$nimages = $e;

# merge them in pairs
$op = "> foo.merge.cat";
for ($l = 0; $l < $nimages; $l++) {
	for ($m = $l + 1; $m < $nimages; $m++) {
		echosys("mergecats $mergetol foo.$l.cat foo.$m.cat $op");
		$op = "| lc -b -o >> foo.merge.cat";
	}
}

echosys("imcattools/mosaics/fitpolymodelmagshift $nimages $fitorder foo -outputarray < foo.merge.cat > foo.A.out");
echosys("cat foo.A.out");

echosys("generatelmodelimage -$halfimsize $halfimsize 256 -$halfimsize $halfimsize 256 < foo.par | contour -g -n 0 -W");

sub echosys {
	warn "$_[0]\n";
	system($_[0]) && die "system($_[0]) failed!\n";
}

