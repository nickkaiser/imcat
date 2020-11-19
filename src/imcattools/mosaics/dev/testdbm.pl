#! /usr/bin/perl

dbmopen(%ASSOC, 'assoc', 0666) || die "Can't open data base!\n";
# $ASSOC{'a'} = "fred";
# $ASSOC{'b'} = "joe";

foreach (keys %ASSOC) {
	print $_, "\t", $ASSOC{$_}, "\n";
}
