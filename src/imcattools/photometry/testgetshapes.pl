#! /usr/bin/perl

# image size
$N=1024;

# psf properties
$a=2.0;
$b=3.0;
$phi=0.0;

# galaxy properties
$ng=100;
$rg=10;
$fg=1000;

# star properties
$ns=100;
$rs=0.001;
$fs=100000;

# background noise level
$sig=1.0;

# detection threshold
$numin = 10.0;

# psfimage size
$M=32;

# applied shear
$gamma = 0.1;



$Mc = $M / 2;
$N2 = 2 * $N;
$a2 = 2 * $a;
$b2 = 2 * $b;


# echosys("ic -c $M $M 'xp $Mc == yp $Mc == mult' | smooth -g $a $b $phi > tmp.psf");
# echosys("make_image 1 $N2 $N2 -f -o $ng $rg $fg -o $ns $rs $fs | transformimage -p 1 $gamma $gamma 1 -c | smooth -g $a2 $b2 $phi ".
#	"| scrunch | ic '%1 grand $sig mult +' - > tmp.fits");
# echosys("hfindpeaks tmp.fits -n $numin | apphot -z 30 | getshapes > tmp.cat");
echosys("plotcat rh mag -x 0 10 25 15 -S < tmp.cat | plotcat e -x -0.5 0.5 -0.5 0.5 -a 1 -S > tmp.stars");
echosys("makestamps -n -M 20.0 tmp.stars");
echosys("combinestamps tmp.stars_stamps/* | spinflip 0 1 -1 0 | cycleimage 1 0 > tmp.psf");



sub echosys {
        print @_, "\n";
        system(@_) && die "$0: System call $_[0] failed!\n";
}

