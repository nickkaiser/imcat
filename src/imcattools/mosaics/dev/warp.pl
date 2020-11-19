#! /usr/bin/perl

# apply quadradic warp

$usage = "warp alpha xcol ycol \napply transformation (x,y) => (x,y) (1 + alpha r^2)\n";

die $usage if ($#ARGV < 0);
$alpha = shift @ARGV;
$xcol = shift @ARGV;
$ycol = shift @ARGV;
while (<>) {
	if (/^#/) {
		print;
	} else {
		chop;
		@F = split(' ');
		$x = $F[$xcol - 1];
		$y = $F[$ycol - 1];
		$rr = $x * $x + $y * $y;
		$F[$xcol - 1] = $x * (1 + $alpha * $rr);
		$F[$ycol - 1] = $y * (1 + $alpha * $rr);
		foreach (@F) {
			printf "%13g", $_;
		}
		print "\n";
	}
}
