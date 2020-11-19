/*
 * iis.h
 *
 * define constants for iis header
 * according to fax from george miyashiro
 */

#define TRANSFER_ID	0
#define THING_COUNT	1
#define	SUB_UNIT	2
#define	CHECK_SUM	3
#define	X_REGISTER	4
#define Y_REGISTER	5
#define Z_REGISTER	6
#define T_REGISTER	7

/* transfer ID definitions */
#define  IREAD             0100000
#define  IWRITE                 00
#define  PACKED             040000
#define  BYPASSIFM          020000
#define  BYTE               010000
#define  ADDWRITE            04000
#define  ACCUM               02000
#define  BLOCKXFER           01000
#define  VRETRACE             0400
#define  MUX32                0200
#define  IMT800               0100

/* Subunits */
#define  REFRESH                 01
#define  LUT                     02
#define  OFM                     03
#define  IFM                     04
#define  FEEDBACK                05
#define  SCROLL                  06
#define  VIDEOM                  07
#define  SUMPROC                 010
#define  GRAPHICS                011
#define  CURSOR                  012
#define  ALU                     013
#define  ZOOM                    014
#define  IPB                     017
/* following from iraf "iis.h" */
#define  IMCURSOR                015
#define  WCS	                 021

/* checksum */
#define	CHECKSUMVAL		0177777

/* from iraf iis.h:
# Command definitions
define  COMMAND           100000B
*/

/* X-register */
#define  ADVXONTC          0100000
#define  ADVXONYOV          040000

/* Y-register */
#define  ADVYONXOV         0100000
#define  ADVYONTC           040000

/* from iraf iis.h:
define  ERASE             100000B               # Erase

# 4 - Button Trackball
define  PUSH               40000B
define  BUTTONA              400B
define  BUTTONB             1000B
define  BUTTONC             2000B
define  BUTTOND             4000B
*/

/* Z-register */
#define  CHAN1                  01
#define  CHAN2                  02
#define  CHAN3                  04
#define  CHAN4                 010
#define  GRCHAN            0100000

/* T-register */
#define  BITPL0                 01
#define  BITPL1                 02
#define  BITPL2                 04
#define  BITPL3                010
#define  BITPL4                020
#define  BITPL5                040
#define  BITPL6               0100
#define  BITPL7               0200
#define  ALLBITPL             0377


/* rest is from iraf iis.h
define  LEN_IISFRAMES           4
define  IISFRAMES       CHAN1, CHAN2, CHAN3, CHAN4

# Colors

define  BLUE                   1B
define  GREEN                  2B
define  RED                    4B
define  MONO                   7B


# IIS Sizes
define  IIS_XDIM              512
define  IIS_YDIM              512
define  MCXSCALE               64       # metacode x scale
define  MCYSCALE               64       # metacode y scale
define  SZB_IISHDR             16       # size of IIS header in bytes
define  SZB_IMCURVAL          160       # size of imcursor value buffer, bytes
define  LEN_ZOOM                3       # zoom parameters
define  LEN_CURSOR              3       # cursor parameters
define  LEN_SPLIT              12       # split screen
define  LEN_LUT               256       # look up table
define  LEN_OFM              1024       # output function look up table
define	SZ_WCSTEXT	      320	# max WCS text chars

# IIS Status Words
define  IIS_FILSIZE             (IIS_XDIM * IIS_YDIM * SZB_CHAR)
define  IIS_BLKSIZE             1024
define  IIS_OPTBUFSIZE          (IIS_XDIM * SZB_CHAR)
define  IIS_MAXBUFSIZE          32768
*/