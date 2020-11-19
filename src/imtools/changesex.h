/* following macro cribbed from DEC's sex.h */
#define swap_half(a) ( ((a & 0xff) << 8) | ((unsigned short)(a) >> 8) )
#define swap_word(a) ( ((a) << 24) | \
                      (((a) << 8) & 0x00ff0000) | \
                      (((a) >> 8) & 0x0000ff00) | \
        ((unsigned int)(a) >>24) )

