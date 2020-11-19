#! /usr/bin/perl

# print albums

require "fields.db";

# $root = "fits";
# $root = "sky";
$root = "sub";

foreach $f (keys %fieldname) {
	$expname = $fieldname{$f};
	&echosys("print_image albums/$expname.$root -f -30 200 > albums/$expname.$root.ps");
	&echosys("lpr -h -s -Pnetps1 albums/$expname.$root.ps");
}
 	

sub echosys {
        print @_, "\n";
        system(@_);
}
