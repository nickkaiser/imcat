#!/usr/bin/perl

# object props
$nobj = 0;
$robj = 3;
$fobj = 1000;

# image size
$N = 256;
$N4 = $N / 4;

# cosmic ray props
$ncr = 3;
$fcr = 1000;

echosys("make_image 1 $N $N -o $nobj $robj $fobj > tmp_obj.fits");

$sigma = 30;
$nim = 5;

for ($i = 0; $i < $nim; $i++) {
	$thesigma = $sigma * ($nim + $i) / $nim;
	echosys("makerandcat $ncr -seed $i | lc x 'f = $fcr' ".
		"| makedensity x 0 1 $N4 0 1 $N4 -v f | unscrunch | unscrunch ". 
		"| ic -s $i -h SIGMA $thesigma '%1 grand $thesigma * + %2 +' tmp_obj.fits - > tmp.$i.fits");
}

sub echosys {
	warn "$_[0]\n";
	system($_[0]) && die "system($_[0]) failed!\n";
}

