/*
 * cat.h
 */


/* types */

#define MAX_DIMS		7

#define	NUM_TYPE		0
#define	TEXT_TYPE		1
#define	COMMENT_TYPE		2
#define	CONTENTS_TYPE		3
#define NUM1_OP_TYPE		4
#define NUM2_OP_TYPE		5
#define GEN_OP_TYPE		6

#define BINARY_FILE_TYPE	0
#define TEXT_FILE_TYPE		1

#define BIG_ENDIAN_CAT_BYTE_ORDER		0
#define LITTLE_ENDIAN_CAT_BYTE_ORDER		1

#ifdef LITTLEENDIAN
#define NATIVE_CAT_BYTE_ORDER		LITTLE_ENDIAN_CAT_BYTE_ORDER
#define NON_NATIVE_CAT_BYTE_ORDER	BIG_ENDIAN_CAT_BYTE_ORDER
#else
#define NATIVE_CAT_BYTE_ORDER		BIG_ENDIAN_CAT_BYTE_ORDER
#define NON_NATIVE_CAT_BYTE_ORDER	LITTLE_ENDIAN_CAT_BYTE_ORDER
#endif


#define	LC_NUM_FMT			" %14.8lg"
#define	LC_TEXT_FMT		" %14s"
#define	LC_RIGHT_TEXT_FMT		" %14s"


#define LC_HEADER_NAME_FMT		"# %-20s %-20s"
#define LC_OBJECT_NAME_FMT		" %s %s"


/* structure definitions */

typedef struct item {
	char 		*name;
	int		itype;
	int		ndim;
	int		idim[MAX_DIMS];
	int		dim[MAX_DIMS];
	void		*addr;
	struct item 	*next;
	int		swapbytesi;
} item;


typedef struct comment {
	char		*text;
	struct comment	*next;
} comment;


typedef	struct cathead {
	item		*headeritembase;
	comment		*commentbase;
	int		nobjectitems;
	item		*objectitembase;
	item		**itemlist;
} cathead;


typedef struct object {
	cathead		*cathead;
	int		nitems;
	void		**addrlist;
	struct object	*next;
} object;



/* function declarations */

cathead		*readcathead(void);
void		copyheaderinfo(cathead *thecathead, cathead *sourcecathead);
void		copycontentinfo(cathead *thecathead, cathead *sourcecathead);
void		writecathead(cathead *thecathead);
void		freecathead(cathead *thecathead);
void		makeitemlist(cathead *thecathead);
int		gettype(void);
item		*readnewitem(int itype);
item		*sscannewitem(int itype, char *string);
item		*newitem(char *name, int itype, int ndim, ...);
item		*newitembydimarray(char *name, int itype, int ndim, int *dim);
item		*copyitem(item *sourceitem);
int		allocitemcontents(item *theitem, void **addr, int level);
void		copyitemcontents(item *dstitem, void *dstaddr, item *srcitem, void *srcaddr, int level);
void		freeitem(item *theitem);
void		freeitemcontents(item *theitem, void *addr, int level);
int		readitem(item *theitem, void *addr, int level);
int		sscanitem(item *theitem, void *addr, int level, char **string);
int		writeitem(item *theitem, void *addr, int level);
int		writelabel(item *theitem, int level);
int		writename(item *theitem);
int		size(int itype);
int		installitem(item *theitem, item **base);
void		installcomment(comment *thecomment, comment **base);
void		addcomment(char *thetext, cathead *thecathead);
void		addargscomment(int argc, char *argv[], cathead *thecathead);
comment		*copycomment(comment *sourcecomment);
object		*newobject(cathead *thecathead);
void		allocobjectcontents(object *theobject);
void		connectcatheadtoobject(object *theobject);
void		connectobjecttocathead(object *theobject);
int		readobject(object *theobject);
void		writeobject(object *theobject);
void		freeobject(object *theobject);
item		*getheaderitem(char *name, cathead *thecathead);
void		*getheaderitemaddress(char *name, cathead *thecathead);
void		setheaderitemaddress(char *name, cathead *thecathead, void *theaddr);
void		deleteheaderitem(char *name, cathead *thecathead);
void		addobjectitem(item *newitem, cathead *thecathead);
void		deleteobjectitem(char *name, cathead *thecathead);
item		*getobjectitem(char *name, cathead *thecathead);
int		getobjectitemindex(char *name, object *theobject);
int		inheritcontents(object *theobject, object *sourceobject);
void		*getaddress(object *theobject, int theindex);
void		setaddress(object *theobject, int theindex, void *theaddr);
void		copystring(char **thestring, char *sourcestring);
void		setcatipf(FILE *stream);
void		setcatopf(FILE *stream);
void		getcatipfiletype(int *filetype);
void		setcatipfiletype(int filetype);
void		getcatopfiletype(int *filetype);
void		setcatopfiletype(int filetype);
void		swapbytes(void *addr, int ndoubles);
void		swapint(void *addr, int nint);
int		install1numericheaderitem(char *name, double value, cathead *thecat);
int		install1textheaderitem(char *name, char *value, cathead *thecat);

/*
COMMENTS:
The general idea is that we want to be able to read in a cat header
which contains various self-defining header items, comments and
a list of self-defining object items and then read a sequence
of objects (one per line) into a suitable object structure.

The types of items are multi-dimensional arrays of types
NUM add TEXT; declared as e.g. "number 2 3 4 F" for
a 2-d array of doubles F[i][j] with i = 0..2, j = 0..3.
Text 'values' are (char *) pointers.

The name, type definition for an item (header item or object item)
is stored in an "item" structure.  This contains name, type, dimensions.
Also contains "void *addr" slot where we can store the value (if it's
a header item - addr slot is left empty for object items).  Also contain
a "item *next" for ease of linking into linked lists.

The header and object item definitions are stored in a "cathead"
structure which contains null-terminated linked lists of header and
object items (so we can easily add, remove items).  It also contains
a "item **itemlist" entry which contains an array of pointers
to the object item addresses (so we can easily retrieve info if we
know the index of the object-item).

We define a "object" structure which contains a pointer to it's
cathead and an array of pointers to the addresses of its contents.

Have also added 2 extra classes of items:

Operators: NUM2_OP_TYPE, NUM1_OP_TYPE, GEN1_OP_TYPE...
These contain
	index to the function (iop) in (theitem->dim)[0]
	address of function in theitem->addr
	space for result in theitem->next;
These are used as stack entries (along with regular items) to
perform calculations.

DESCRIPTION OF FUNCTIONS:

cathead		*readcathead(void);
	Creates a new cathead and reads it from the catipf stream.

void		copyheaderinfo(cathead *cathead, cathead *sourcecathead);
	Fill empty cathead with duplicates of sourcecathead's header items
	and comments.

void		copycontentinfo(cathead *cathead, cathead *sourcecathead);
	Make a duplicate of linked list headed by objectitembase and
	also its addrlist.

void		writecathead(cathead *thecathead);
	Writes a cathead to catopf.

void		freecathead(cathead *thecathead);
	Frees a cathead and all of it's contents.

void		makeitemlist(cathead *thecathead);
	Construct the list of object items (trashes the old one if
	it already exists).  Called by readcathead(), addobjectitem(),
	deleteobjectitem().

int		gettype(void);
	Reads a word from catipf (should be int, float etc) and returns
	the integer "itype".

item		*readnewitem(int itype);
	Creates a new item by reading ndim, dim[0], dim[1],.... from catipf.

item		*newitem(char *name, int itype, int ndim, ...);
item		*newitembydimarray(char *name, int itype, int ndim, int *dim);
	Creates a new item of specified name, type, size with either varargs method
	or by passing dimensions as an array.

item		*copyitem(item *sourceitem);
	Creates a copy of "sourceitem"
	Also allocates new space for "name" if it exists, but does
	try to copy the contents (item->addr and item->next are set to NULL.

int		allocitemcontents(item *theitem, void **addr, int level);
	Allocates memory for contents of an item under "addr", using
	template in "theitem". Works recursively.

void		freeitem(item *theitem);
	Frees up an item.  Calls freeitemcontents().

void		freeitemcontents(item *theitem, void *addr, int level);
	Free's contents stored in "addr". Works recursively.

int		readitem(item *theitem, void *addr, int level);
	Read contents of an item from catipf stream recursively.

int		writeitem(item *theitem, void *addr, int level);
	Write contents of an item to catopf stream recursively.

int		writelabel(item *theitem, int level);
	Writes a lable for an item.  Called by writecathead().

int		writename(item *theitem);
	Writes type, sizeinfo and name of an item to catopf.

int		size(int itype);
	Get the size required to store a given type.

int		installitem(item *theitem, item **base);
	Install an item at the end of a null-terminated linked list.
	Returns 1 if it replaced an existing item, 0 otherwise.

void		installcomment(comment *thecomment, comment **base);
	Install a comment at the end of a null-terminated linked list.

void		addcomment(char *thetext, cathead *thecathead);
void		addargscomment(int argc, char *argv[], cathead *thecathead);
	Add comments.

comment		*copycomment(comment *sourcecomment);
	Generate a new comment and duplicate contents of sourcecomment.
	Makes new space for the "text".

object		*newobject(cathead *thecathead);
	Creates a new object using info in cathead.

void		allocobjectcontents(object *theobject);
	Allocate space for any empty addr[index].

void		connectcatheadtoobject(object *theobject);
	Makes the cathead itemlist addresses point at the object.

void		connectobjecttocathead(object *theobject);
	Make the object addrlist entries point at the items in the cathead.

int		readobject(object *theobject);
void		writeobject(object *theobject);
	Read/write contents of an obect from/to catipf/catopf.

void		freeobject(object *theobject);
	Free an object and all of its contents.

item		*getheaderitem(char *name, cathead *thecathead);
void		*getheaderitemaddress(char *name, cathead *thecathead);
void		setheaderitemaddress(char *name, cathead *thecathead, void *theaddr);
	Get/set addresses for contents of header items.

void		deleteheaderitem(char *name, cathead *thecathead);
	Remove a header item.

void		addobjectitem(item *newitem, cathead *thecathead);
void		deleteobjectitem(char *name, cathead *thecathead);
	Add/delete object items from cathead.  Calls makeitemlist();

item		*getobjectitem(char *name, cathead *thecathead);
	Get object item by name.

int		getobjectitemindex(char *name, object *theobject);
	Get the index for an object item by name so we can dereference
	efficiently.

int		inheritcontents(object *theobject, object *sourceobject);
	Copy all the addresses from "sourceobject" which are also named
	items in "theobject".  Any other items in theobject will be empty.

void		*getaddress(object *theobject, int theindex);
void		setaddress(object *theobject, int theindex, void *theaddr);
	Get/set addresses for contents of object by index.

void		copystring(char **thestring, char *sourcestring);
	Make a duplicate of a string.  Free's "thestring" if non-empty.

void		setcatipf(FILE *stream);
void		setcatopf(FILE *stream);
	Change the file input/output streams.

void		getcatipfiletype(int *filetype);
void		setcatipfiletype(int filetype);
void		getcatopfiletype(int *filetype);
void		setcatopfiletype(int filetype);
	get and set the input and output file types

void		swapbytes(void *addr, int ndoubles);
void		swapint(void *addr, int nint);
	for reading/writing in non-native format

int		install1numericheaderitem(char *name, double value, cathead *thecat);
int		install1textheaderitem(char *name, char *value, cathead *thecat);
	for installing single header items
*/
