#! /usr/bin/perl

# to generate convertarray.c --- convert 1-D array from e.g. float to int

$ntypes = 5;
@pixtype = ("UCHAR_PIXTYPE", "SHORT_PIXTYPE", "INT_PIXTYPE", "FLOAT_PIXTYPE", "DBL_PIXTYPE");
@type =    ("unsigned char", "short", "int", "float", "double");
@conv =    ("(unsigned char) clip(0, UCHAR_MAX - 1, floor(0.5 +",
		"(short) clip(SHRT_MIN + 1, SHRT_MAX, floor(0.5 +", 
		"(int) clip(INT_MIN + 1, INT_MAX, floor(0.5 +", 
		"(float) ((", 
		"(double) ((");

@magic =   ("UCHAR_MAGIC", "SHORT_MAGIC", "INT_MAGIC", "FLOAT_MAGIC", "DBL_MAGIC");


$declaration = "int\tconvertarray(char *fsrc, char *fdst, int srcpixtype, int dstpixtype, int nel, int bscaling, double bscale, double bzero)";

open(OPF, ">convertarray.h");
print OPF "/*\n * convertarray.h -- written by makeconvertarray.pl\n */\n\n$declaration;\n";
close(OPF);

open(OPF, ">convertarray.c");

print OPF "/*\n * convertarray.c -- written by makeconvertarray.pl\n */\n\n";

print OPF "#include <stdio.h>\n#include <math.h>\n#include <limits.h>\n#include <values.h>\n";
print OPF "#include \"fits.h\"\n#include \"convertarray.h\"\n\n";

print OPF "#define clip(min,max,f) ((f) > (min) ? ((f) < (max) ? (f) : (max)) : (min))\n\n";

print OPF "$declaration\n{\n";
print OPF "\tint\ti, srcsize, dstsize;\n\n";

print OPF "\tswitch(srcpixtype) {\n";
for ($t = 0; $t < $ntypes; $t++) {
	print OPF "\t\tcase $pixtype[$t]:\n";
	print OPF "\t\t\tsrcsize = sizeof($type[$t]);\n";
	print OPF "\t\t\tbreak;\n";
}
print OPF "\t\tdefault:\n";
print OPF "\t\t\terror_exit(\"convertarray: bad pixtype\\n\");\n";
print OPF "\t}\n\n";

print OPF "\tswitch(dstpixtype) {\n";
for ($t = 0; $t < $ntypes; $t++) {
	print OPF "\t\tcase $pixtype[$t]:\n";
	print OPF "\t\t\tdstsize = sizeof($type[$t]);\n";
	print OPF "\t\t\tswitch(srcpixtype) {\n";
	for ($tt = 0; $tt < $ntypes; $tt++) {
		print OPF "\t\t\t\tcase $pixtype[$tt]:\n";
		print OPF "\t\t\t\t\tfor (i = 0; i < nel; i++) {\n";
		print OPF "\t\t\t\t\t\tif (*(($type[$tt] *) (fsrc + i * srcsize)) == $magic[$tt]) {\n";
		print OPF "\t\t\t\t\t\t\t*(($type[$t] *) (fdst + i * dstsize)) = $magic[$t];\n";
		print OPF "\t\t\t\t\t\t} else {\n";
		print OPF "\t\t\t\t\t\t\tif (bscaling) {\n";
		print OPF "\t\t\t\t\t\t\t\t*(($type[$t] *) (fdst + i * dstsize)) = $conv[$t] *(($type[$tt] *) (fsrc + i * srcsize)) * bscale + bzero));\n";
		print OPF "\t\t\t\t\t\t\t} else {\n";
		print OPF "\t\t\t\t\t\t\t\t*(($type[$t] *) (fdst + i * dstsize)) = $conv[$t] *(($type[$tt] *) (fsrc + i * srcsize))));\n";
		print OPF "\t\t\t\t\t\t\t}\n";
		print OPF "\t\t\t\t\t\t}\n";
		print OPF "\t\t\t\t\t}\n";
		print OPF "\t\t\t\t\tbreak;\n";
	}
	print OPF "\t\t\t\tdefault:\n";
	print OPF "\t\t\t\t\terror_exit(\"convertarray: bad pixtype\\n\");\n";
	print OPF "\t\t\t}\n";
	print OPF "\t\t\tbreak;\n";
}
print OPF "\t\tdefault:\n";
print OPF "\t\t\terror_exit(\"convertarray: bad pixtype\\n\");\n";
print OPF "\t}\n\n";



print OPF "}\n";
close(OPF);
