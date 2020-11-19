
#define mpegparamtext "# parameter file template with lots of comments to assist you\n\
#\n\
# you can use this as a template, copying it to a separate file then modifying\n\
# the copy\n\
#\n\
#\n\
# any line beginning with '#' is a comment\n\
#\n\
# no line should be longer than 255 characters\n\
#\n\
#\n\
# general format of each line is:\n\
#	<option> <spaces and/or tabs> <value>\n\
#\n\
# lines can generally be in any order\n\
#\n\
# an exception is the option 'INPUT' which must be followed by input\n\
# files in the order in which they must appear, followed by 'END_INPUT'\n\
#\n\
# Also, if you use the `command` method of generating input file names,\n\
# the command will only be executed in the INPUT_DIR if INPUT_DIR preceeds\n\
# the INPUT parameter.\n\
#\n\
# <option> MUST be in UPPER CASE\n\
#\n\
\n\
PATTERN		IBBPBBPBBPBBPBBP\n\
OUTPUT		tmp.mpg\n\
\n\
# mpeg_encode really only accepts 3 different file formats, but using a\n\
# conversion statement it can effectively handle ANY file format\n\
#\n\
# You must specify the type of the input files.  The choices are:\n\
#    YUV, PPM, JMOVIE, Y, JPEG, PNM\n\
#	(must be upper case)\n\
#\n\
BASE_FILE_FORMAT	PPM\n\
\n\
#\n\
# if YUV format (or using parallel version), must provide width and height\n\
# YUV_SIZE	widthxheight\n\
# this option is ignored if BASE_FILE_FORMAT is not YUV and you're running\n\
# on just one machine\n\
#\n\
YUV_SIZE	352x240\n\
\n\
# If you are using YUV, there are different supported file formats.\n\
# EYUV or UCB are the same as previous versions of this encoder.\n\
# (All the Y's, then U's then V's, in 4:2:0 subsampling.)\n\
# Other formats, such as Abekas, Phillips, or a general format are\n\
# permissible, the general format is a string of Y's, U's, and V's\n\
# to specify the file order.\n\
\n\
INPUT_FORMAT UCB\n\
\n\
# the conversion statement\n\
#\n\
# Each occurrence of '*' will be replaced by the input file\n\
#\n\
# e.g., if you have a bunch of GIF files, then this might be:\n\
#	INPUT_CONVERT	giftoppm *\n\
#\n\
# e.g., if you have a bunch of files like a.Y a.U a.V, etc., then:\n\
#	INPUT_CONVERT	cat *.Y *.U *.V\n\
#\n\
# e.g., if you are grabbing from laser disc you might have something like\n\
#	INPUT_CONVERT	goto frame *; grabppm\n\
# 'INPUT_CONVERT *' means the files are already in the base file format\n\
#\n\
INPUT_CONVERT	*\n\
\n\
# number of frames in a GOP.\n\
#\n\
# since each GOP must have at least one I-frame, the encoder will find the\n\
# the first I-frame after GOP_SIZE frames to start the next GOP\n\
#\n\
# later, will add more flexible GOP signalling\n\
#\n\
GOP_SIZE	16\n\
\n\
# number of slices in a frame\n\
#\n\
# 1 is a good number.  another possibility is the number of macroblock rows\n\
# (which is the height divided by 16)\n\
#\n\
SLICES_PER_FRAME	1\n\
\n\
# directory to get all input files from (makes this file easier to read)\n\
INPUT_DIR	./tmp\n\
\n\
# There are a bunch of ways to specify the input files.\n\
# from a simple one-per-line listing, to the following \n\
# way of numbering them.  See the manual for more information.\n\
INPUT\n\
# '*' is replaced by the numbers 01, 02, 03, 04\n\
# if I instead do [01-11], it would be 01, 02, ..., 09, 10, 11\n\
# if I instead do [1-11], it would be 1, 2, 3, ..., 9, 10, 11\n\
# if I instead do [1-11+3], it would be 1, 4, 7, 10\n\
# the program assumes none of your input files has a name ending in ']'\n\
# if you do, too bad!!!\n\
#\n\
#\n\
*.ppm	[00000-99999]\n\
# can have more files here if you want...there is no limit on the number\n\
# of files\n\
END_INPUT\n\
\n\
\n\
\n\
# Many of the remaining options have to do with the motion search and qscale\n\
\n\
# FULL or HALF -- must be upper case\n\
PIXEL		HALF\n\
\n\
# means +/- this many pixels for both P and B frame searches\n\
# specify two numbers if you wish to serc different ranges in the two.\n\
RANGE		10\n\
\n\
# this must be one of {EXHAUSTIVE, SUBSAMPLE, LOGARITHMIC}\n\
PSEARCH_ALG	LOGARITHMIC\n\
\n\
# this must be one of {SIMPLE, CROSS2, EXHAUSTIVE}\n\
#\n\
# note that EXHAUSTIVE is really, really, really slow\n\
#\n\
BSEARCH_ALG	CROSS2\n\
\n\
#\n\
# these specify the q-scale for I, P, and B frames\n\
# (values must be between 1 and 31)\n\
# These are the Qscale values for the entire frame in variable bit-rate\n\
# mode, and starting points (but not important) for constant bit rate\n\
#\n\
IQSCALE		8\n\
PQSCALE		10\n\
BQSCALE		25\n\
\n\
# this must be ORIGINAL or DECODED\n\
REFERENCE_FRAME	ORIGINAL\n\
\n\
# for parallel parameters see parallel.param in the exmaples subdirectory\n\
\n\
# if you want constant bit-rate mode, specify it as follows (number is bits/sec):\n\
BIT_RATE  1000000\n\
\n\
# To specify the buffer size (327680 is default, measused in bits, for 16bit words)\n\
BUFFER_SIZE 327680\n\
\n\
# The frame rate is the number of frames/second (legal values:\n\
# 23.976, 24, 25, 29.97, 30, 50 ,59.94, 60\n\
FRAME_RATE 30\n\
\n\
# There are many more options, see the users manual for examples....\n\
# ASPECT_RATIO, USER_DATA, GAMMA, IQTABLE, etc.\n"

