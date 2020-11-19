#ifdef LITTLEENDIAN
        #define _D0     3
#else
        #define _D0     0
#endif

#define _DOFF   4
#define _DMAX   ((1<<(15-_DOFF))-1)
#define _DNAN   (0x8000|_DMAX<<_DOFF|1<<(_DOFF-1))

#define NBITS   (48+_DOFF)
#if _D0
#define INIT(w0)        0, 0, 0, w0
#else
#define INIT(w0)        w0, 0, 0, 0
#endif

typedef const union {
        unsigned short _W[4];
        double  _D;
} _Dconst;

