#!/usr/bin/perl

# now we have run mosaicfit to get 1st approcimation to mosaicfit.par
# we want to make a refined mergelist2.out

require "fields.db";
$sldir = "superlists";

&echosys("rm mergelists2.out");
foreach $f (keys %fieldname) {
	&echosys("transformsuperlist.pl < $sldir/field$f.sl > temp$f.sl");
}
foreach $f1 (keys %fieldname) {
	foreach $f2 (keys %fieldname) {
		next if ($f2 <= $f1);
		&echosys("mergelists temp$f1.sl temp$f2.sl -d 5 > temp.merge");
		&echosys("cat  temp.merge >> mergelists2.out");
	}
}
&echosys("doto mergelists2.out transformmergedlist.pl");
&echosys("rm temp.sl temp.merge");


sub echosys {
        print @_, "\n";
        system(@_);
}
