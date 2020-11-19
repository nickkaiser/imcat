void		dosky(				float **f, 
								int N1, 
								int N2, 
								int i, 
								int j, 
								int a1, 
								int a2, 
								skyquad * sky, 
								float fmode,
								float sigma);
void		do_object_stats(	object *pk, 
					float **f, 
					int N1, 
					int N2, 
					float(*fsky)(int i, int j), 
					float sigma,
					float ne,
					float rc,
					float alpha);
int			objectcmp(object *pk1, object *pk2);
int			shortcmp(short *s1, short *s2);
float		fsky(int di, int dj);
void		setskyparameters(skyquad *sky);
void		zap(int zapmode, object *obj, int radiustype, float a, float **f, float **fzap, short **nzap, int N1, int N2);






