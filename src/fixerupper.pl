#!/usr/bin/perl

system("grep usr/bin/perl `find . -type f` > tmp.lst");
@scripts = split("\n", `grep -v README < tmp.lst | grep -v Binary`);
foreach $script (@scripts) {
	($script) = split(":", $script);
	print "$script\n";
	$com = "mv $script tmp.pl ; sed 's:usr:usr:' < tmp.pl > $script ; chmod +x $script";
	print "$com\n";
	system($com);
}


exit;
