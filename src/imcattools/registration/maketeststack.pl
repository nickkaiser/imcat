#!/usr/bin/perl

$| = 1;

$Ncats = 9;
$Nobjects = 20;
$dmax = 100;
$dphimax = 0.02;
$alpha = 3.e-8;
$dmagmax = 1.0;

$catdir = "images";
$basecat = "$catdir/base.cat";

open(EXPDB, ">exp.db");

open(PARFILE, ">maketeststack.par");
printf PARFILE "#          dx            dy         phi00         phi01         phi10         phi11            dm\n";

&echosys("maketestcat $Nobjects 3 | lc -C -N '1 2 x' -n mag | lc +all 'x = %x 2048 vscale' > $basecat");
for ($i = 0; $i < $Ncats; $i++) {
	$dstcat = "$catdir/".$i.".cat";
	if ($i) {
		$dx = 2 * $dmax * (rand() - 0.5);
		$dy = 2 * $dmax * (rand() - 0.5);
		$phi00 = 1.0 + 2 * $dphimax * (rand() - 0.5);
		$phi11 = 1.0 + 2 * $dphimax * (rand() - 0.5);
		$phi01 = 2 * $dphimax * (rand() - 0.5);
		$phi10 = 2 * $dphimax * (rand() - 0.5);
		$dmag = 2 * $dmagmax * (rand() - 0.5);
	} else {
		$phi00 = $phi11 = 1.0;
		$dx = $dy = $phi01 = $phi10 = $dmag = 0.0;
	}
	printf PARFILE "%13f %13f %13f %13f %13f %13f %13f\n", $dx, $dy, 2 - $phi00, -$phi01, -$phi10, 2 - $phi11, $dmag;
	&echosys("lc +all 'x = %x $dx $dy 2 vector vsub $phi00 $phi01 $phi10 $phi11 lintrans' 'mag = %mag $dmag +' 'sym = $i' < $basecat > $dstcat");
	$srccat = $dstcat;
	$dstcat = "$catdir/dist_".$i.".cat";
	&echosys("lc +all 'x = %x -1024.0 vshift enter enter dot $alpha mult 1 + vscale 1024.0 vshift' < $srccat > $dstcat");
#	print EXPDB "$i\n";
	print EXPDB "dist_$i\n";
}
close(EXPDB);
close(PARFILE);



sub echosys {
        print @_, "\n";
        system(@_) && die "System call failed!\n";
}
