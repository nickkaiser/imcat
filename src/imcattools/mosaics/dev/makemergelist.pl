#!/usr/bin/perl


require "fields.db";
$sldir = "superlists";

&echosys("rm mergelists.out");
foreach $exp1 (keys %fieldname) {
	foreach $exp2 (keys %fieldname) {
		next if ($exp2 <= $exp1);
		print $exp1, $exp2, "\n";
		print "# calculating transformation....\n";
		&echosys("register.scr $sldir/field$exp1.sl $sldir/field$exp2.sl  256 30 > temp.par");
		system("cat temp.par");
		&echosys("scalerottrans `cat temp.par` < $sldir/field$exp2.sl > temp.sl");
		&echosys("mergelists $sldir/field$exp1.sl temp.sl -d 30 > temp.merge");
		&echosys("scalerottrans `cat temp.par` -i -c 5 6 < temp.merge >> mergelists.out");
	}
}
&echosys("rm temp.sl temp.merge");


sub echosys {
        print @_, "\n";
        system(@_);
}
