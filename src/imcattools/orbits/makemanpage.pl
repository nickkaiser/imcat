#!/usr/bin/perl

$usage = "usage : $0 command [-print]\n";

$tool = shift(@ARGV) || die $usage;
$print = 1 if ($#ARGV >= 0);

echosys("$tool -u 2> tmp.man");
if ($print) {
	echosys("lpr tmp.man ; rm tmp.man");
} else {
	echosys("makegroffmanpage $tool 'Celestial_Mechanics' < tmp.man > ~/man/man1/$tool.1 ; rm tmp.man");
}

sub echosys {
        warn "$_[0]\n";
        system($_[0]);
}
