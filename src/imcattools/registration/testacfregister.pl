#!/usr/bin/perl

$| = 1;

$a = 1.3;
$phi = 0.23;
$dx = 120.0;
$dy = 260.0;
$N = 2048;
$nobj = 20;

print "# input parameters: $dx $dy $a $phi\n";
echosys("makerandcat $nobj | lc 'x = %x $N vscale' > tmp_0.cat");
echosys("scalerottrans $dx $dy $a $phi < tmp_0.cat > tmp_1.cat");
echosys("acfregister tmp_[01].cat");


sub echosys {
        print @_, "\n";
        system(@_) && die "System call failed!\n";
}
