#! /usr/bin/perl

# take a random sample

$usage = "randomsample p\nchoose objects at random with probability p\n";

die $usage if ($#ARGV < 0);
$p = shift @ARGV;
srand;
while (<>) {
	if (/^#/) {
		print;
	} else {
		print $_ if (rand() < $p);
	}
}
