/* 
 * CDLPROTO.H  -- CDL Library internal function prototype definitions.
 */


/* IMD Interface Function definitions.  */

#ifdef ANSI_FUNC

IMDPtr imd_open(char *imtdev);
int imd_displayImage(IMDPtr imd, uchar *pix, int nx, int ny, int frame, int fbconfig, int comp_wcs);
int imd_readCursor(IMDPtr imd, int sample, float *x, float *y, char *key);
int imd_setWCS(IMDPtr imd, char *name, char *title, float a, float b, float c, float d, float tx, float ty, float z1, float z2, int zt);
int imd_getWCS(IMDPtr imd, char *name, char *title, float *a, float *b, float *c, float *d, float *tx, float *ty, float *z1, float *z2, int *zt);
int imd_close(IMDPtr imd);
int imd_writeImage(IMDPtr imd, uchar *pix, int nx, int ny, int llx, int lly);
int imd_readImage(IMDPtr imd, uchar *pix, int *nx, int *ny);
int imd_readFrameBuffer(IMDPtr imd, uchar *pix, int *nx, int *ny);
int imd_setFrame(IMDPtr imd, int frame);
int imd_setFBConfig(IMDPtr imd, int configno);
int imd_getFBConfig(IMDPtr imd, int *configno, int *width, int *height, int *nframes);
int imd_setName(IMDPtr imd, char *name);
int imd_setTitle(IMDPtr imd, char *title);
int imd_setCursor(IMDPtr imd, int x, int y, int wcs);
int imd_clearFrame(IMDPtr imd);
int imd_readSubRaster(IMDPtr imd, int llx, int lly, int nx, int ny, uchar *pix);
int imd_writeSubRaster(IMDPtr imd, int llx, int lly, int nx, int ny, uchar *pix);
int imd_setDebug(int state);

#else
                     
IMDPtr  	imd_open();
int     	imd_setFBConfig(), imd_writeDisplay(), imd_readDisplay();
int		imd_setFrame(), imd_setCursor(), imd_readCursor();
int     	imd_close(), imd_readSubRaster(), imd_writeSubRaster();
int     	imd_clearFrame(), imd_setWCS(), imd_getWCS();

#endif



/* COMM Interface Function definitions.  */

#ifdef ANSI_FUNC

int com_writeData(int fd, short x, short y, uchar *pix, int nbytes);
int com_readData(int fdin, int fdout, short x, short y, uchar *pix, int *npix);
int com_readCursor(int fdin, int fdout, int sample, float *x, float *y, char *key);
int com_setCursor(int fd, int x, int y, int wcs);
int com_setFBConfig(int fd, int configno);
int com_setFrame(int fd, int frame_num);
int com_writeWCS(int fd, char *name, float a, float b, float c, float d, float tx, float ty, float z1, float z2, int zt);
int com_readWCS(int fdin, int fdout, char *name, float *a, float *b, float *c, float *d, float *tx, float *ty, float *z1, float *z2, int *zt);
int com_eraseFrame(int fd);

int com_setDebug(int state);

#else

int		com_writeData(), com_readData(), com_setFBConfig();
int		com_readCursor(), com_setCursor(), com_setFrame();
int		com_writeWCS(), com_readWCS(), com_eraseFrame();
int 		com_setDebug();

#endif



/* EPS Interface Function definitions.  */

#ifdef ANSI_FUNC

PSImage *eps_init(void);
void eps_print(PSImage *psim, FILE *fp, uchar *data, int xdim, int ydim, int depth, int pad);
void eps_close(PSImage *psim);
void eps_setPage(PSImage *psim, int orientation, int paper_size, int scale, int flags);
void eps_setCmap(PSImage *psim, uchar *r, uchar *g, uchar *b, int ncolors);
void eps_setCompress(PSImage *psim, int compress);
void eps_setColorType(PSImage *psim, int color_class);
void eps_setLabel(register PSImage *psim, char *label);
void eps_setTransform(PSImage *psim, float z1, float z2, int ztype, float offset, float scale, char *cmap_name);
void eps_setCorners(PSImage *psim, int llx, int lly, int urx, int ury);
void eps_getImageSize(PSImagePtr psim, int xdim, int ydim, float *width, float *height);
void eps_getImagePos(PSImagePtr psim, int xdim, int ydim, int *llx, int *lly);

#else

PSImage *eps_init();
void 	eps_print(), eps_close(), eps_setPage(), eps_setCmap();
void 	eps_setCompress(), eps_setColorType(), eps_setLabel();
void	eps_setTransform(), eps_setCorners(), eps_getImageSize();
void	eps_getImagePos();

#endif
