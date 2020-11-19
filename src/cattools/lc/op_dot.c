/*
 * op_dot.c - vector dot product
 */

#include <stdio.h>
#include "../../catlib/cat.h"
#include "stack.h"
#include "error.h"

void	dotinit(item *operator)
{
	item	*result, *op1, *op2;
	int	i, j, k, ndim1, ndim2, dimsum;
	double	*product, *x, *y, **m1, **m2, **mproduct;

	op2 = pop();
	op1 = pop();
	ndim1 = op1->ndim;
	ndim2 = op2->ndim;

	if ((op1->itype != NUM_TYPE) || (op2->itype != NUM_TYPE) ||
		((op1->dim)[op1->ndim - 1] != (op2->dim)[0])) {
		error_exit("dot: non-compatible arguments\n");
	}
	if ((ndim1 == 1) && (ndim2 = 1)) {	/* vector-vector product */
		result = operator->next = newitem("DOT_TEMP", NUM_TYPE, 1, 1);
		product = (double *) calloc(1, sizeof(double));
		result->addr = (void *) product;
		x = (double *) (op1->addr);
		y = (double *) (op2->addr);
		for (i = 0; i < (op1->dim)[0]; i++) {
			*product += x[i] * y[i];
		}
	} else if ((ndim1 == 2) && (ndim2 == 1)) {  /* matrix-vector product */
		result = operator->next = newitem("DOT_TEMP", NUM_TYPE, 1, op1->dim[0]);
		product = (double *) calloc(op1->dim[0], sizeof(double));
		result->addr = (void *) product;
		m1 = (double **) (op1->addr);
		y = (double *) (op2->addr);
		for (j = 0; j < op1->dim[0]; j++) {
			for (i = 0; i < (op2->dim)[0]; i++) {
				product[j] += m1[j][i] * y[i];
			}
		}
	} else if ((ndim1 == 2) && (ndim2 == 2)) {  /* matrix-matrix product */
		result = operator->next = newitem("DOT_TEMP", NUM_TYPE, 2, op1->dim[0], op2->dim[ndim2-1]);
		mproduct = (double **) calloc(op1->dim[0], sizeof(double *));
		for (j = 0; j < op1->dim[0]; j++) {
			mproduct[j] = (double *) calloc(op2->dim[ndim2 - 1], sizeof(double));
		}
		result->addr = (void *) mproduct;
		m1 = (double **) (op1->addr);
		m2 = (double **) (op2->addr);
		for (j = 0; j < op1->dim[0]; j++) {
			for (k = 0; k < op2->dim[ndim2-1]; k++) {
				for (i = 0; i < (op2->dim)[0]; i++) {
					mproduct[j][k] += m1[j][i] * m2[i][k];
				}
			}
		}
	} else {
		error_exit("dot: illegal dimensions\n");
	}
}

void	dotdoit(item *operator)
{
	item	*op1, *op2;
	int	i, j, k, ndim1, ndim2;
	double	*product, *x, *y, **m1, **m2, **mproduct;

	op2 = pop();
	op1 = pop();
	ndim1 = op1->ndim;
	ndim2 = op2->ndim;

	if ((ndim1 == 1) && (ndim2 = 1)) {	/* vector-vector product */
		product = (double *) ((operator->next)->addr);
		x = (double *) (op1->addr);
		y = (double *) (op2->addr);
		*product = 0.0;
		for (i = 0; i < (op1->dim)[0]; i++) {
			*product += x[i] * y[i];
		}
	} else if ((ndim1 == 2) && (ndim2 == 1)) {  /* matrix-vector product */
		product = (double *) ((operator->next)->addr);
		m1 = (double **) (op1->addr);
		y = (double *) (op2->addr);
		for (j = 0; j < op1->dim[0]; j++) {
			product[j] = 0.0;
			for (i = 0; i < (op2->dim)[0]; i++) {
				product[j] += m1[j][i] * y[i];
			}
		}
	} else if ((ndim1 == 2) && (ndim2 == 2)) {  /* matrix-matrix product */
		mproduct = (double **) ((operator->next)->addr);
		m1 = (double **) (op1->addr);
		m2 = (double **) (op2->addr);
		for (j = 0; j < op1->dim[0]; j++) {
			for (k = 0; k < op2->dim[ndim2-1]; k++) {
				mproduct[j][k] = 0.0;
				for (i = 0; i < (op2->dim)[0]; i++) {
					mproduct[j][k] += m1[j][i] * m2[i][k];
				}
			}
		}
	} else {
		error_exit("dot: illegal dimensions\n");
	}

}


void	vsubinit(item *operator)
{
	item	*result, *op1, *op2;
	int	i;
	double	*res, *x, *y;

	op2 = pop();
	op1 = pop();
	if ((op1->itype != NUM_TYPE) || (op2->itype != NUM_TYPE) ||
		(op1->ndim != 1) || (op2->ndim != 1) ||
		((op1->dim)[0] != (op2->dim)[0])) {
		error_exit("dot: non-identical or non-numeric arguments\n");
	}
	result = operator->next = newitem("VSUB_TEMP", NUM_TYPE, 1, (op1->dim)[0]);
	res = (double *) calloc((op1->dim)[0], sizeof(double));
	result->addr = (void *) res;
	x = (double *) (op1->addr);
	y = (double *) (op2->addr);
	for (i = 0; i < (op1->dim)[0]; i++) {
		res[i] = x[i] - y[i];
	}
}

void	vsubdoit(item *operator)
{
	item	*op1, *op2;
	int	i;
	double	*res, *x, *y;

	op2 = pop();
	op1 = pop();
	res = (double *) ((operator->next)->addr);
	x = (double *) (op1->addr);
	y = (double *) (op2->addr);
	for (i = 0; i < (op1->dim)[0]; i++) {
		res[i] = x[i] - y[i];
	}
}

void	msubinit(item *operator)
{
	item	*result, *op1, *op2;
	int	i, j;
	double	**res, **x, **y;

	op2 = pop();
	op1 = pop();
	if ((op1->itype != NUM_TYPE) || (op2->itype != NUM_TYPE) ||
		(op1->ndim != 2) || (op2->ndim != 2) ||
		((op1->dim)[0] != (op2->dim)[0]) ||
		((op1->dim)[1] != (op2->dim)[1])) {
		error_exit("msubinit: non-identical or non-numeric arguments\n");
	}
	result = operator->next = newitem("MSUB_TEMP", NUM_TYPE, 2, (op1->dim)[0], (op1->dim)[1]);
	res = (double **) calloc((op1->dim)[0], sizeof(double *));
	for (i = 0; i < (op1->dim)[0]; i++) {
		res[i] = (double *) calloc((op1->dim)[1], sizeof(double));
	}
	result->addr = (void *) res;
	x = (double **) (op1->addr);
	y = (double **) (op2->addr);
	for (i = 0; i < (op1->dim)[0]; i++) {
		for (j = 0; j < (op1->dim)[1]; j++) {
			res[i][j] = x[i][j] - y[i][j];
		}
	}
}

void	msubdoit(item *operator)
{
	item	*op1, *op2;
	int	i, j;
	double	**res, **x, **y;

	op2 = pop();
	op1 = pop();
	res = (double **) ((operator->next)->addr);
	x = (double **) (op1->addr);
	y = (double **) (op2->addr);
	for (i = 0; i < (op1->dim)[0]; i++) {
		for (j = 0; j < (op1->dim)[1]; j++) {
			res[i][j] = x[i][j] - y[i][j];
		}
	}
}

void	vaddinit(item *operator)
{
	item	*result, *op1, *op2;
	int	i;
	double	*res, *x, *y;

	op2 = pop();
	op1 = pop();
	if ((op1->itype != NUM_TYPE) || (op2->itype != NUM_TYPE) ||
		(op1->ndim != 1) || (op2->ndim != 1) ||
		((op1->dim)[0] != (op2->dim)[0])) {
		error_exit("dot: non-identical or non-numeric arguments\n");
	}
	result = operator->next = newitem("VSUB_TEMP", NUM_TYPE, 1, (op1->dim)[0]);
	res = (double *) calloc((op1->dim)[0], sizeof(double));
	result->addr = (void *) res;
	x = (double *) (op1->addr);
	y = (double *) (op2->addr);
	for (i = 0; i < (op1->dim)[0]; i++) {
		res[i] = x[i] + y[i];
	}
}

void	vadddoit(item *operator)
{
	item	*op1, *op2;
	int	i;
	double	*res, *x, *y;

	op2 = pop();
	op1 = pop();
	res = (double *) ((operator->next)->addr);
	x = (double *) (op1->addr);
	y = (double *) (op2->addr);
	for (i = 0; i < (op1->dim)[0]; i++) {
		res[i] = x[i] + y[i];
	}
}

void	maddinit(item *operator)
{
	item	*result, *op1, *op2;
	int	i, j;
	double	**res, **x, **y;

	op2 = pop();
	op1 = pop();
	if ((op1->itype != NUM_TYPE) || (op2->itype != NUM_TYPE) ||
		(op1->ndim != 2) || (op2->ndim != 2) ||
		((op1->dim)[0] != (op2->dim)[0]) ||
		((op1->dim)[1] != (op2->dim)[1])) {
		error_exit("msubinit: non-identical or non-numeric arguments\n");
	}
	result = operator->next = newitem("MADD_TEMP", NUM_TYPE, 2, (op1->dim)[0], (op1->dim)[1]);
	res = (double **) calloc((op1->dim)[0], sizeof(double *));
	for (i = 0; i < (op1->dim)[0]; i++) {
		res[i] = (double *) calloc((op1->dim)[1], sizeof(double));
	}
	result->addr = (void *) res;
	x = (double **) (op1->addr);
	y = (double **) (op2->addr);
	for (i = 0; i < (op1->dim)[0]; i++) {
		for (j = 0; j < (op1->dim)[1]; j++) {
			res[i][j] = x[i][j] + y[i][j];
		}
	}
}

void	madddoit(item *operator)
{
	item	*op1, *op2;
	int	i, j;
	double	**res, **x, **y;

	op2 = pop();
	op1 = pop();
	res = (double **) ((operator->next)->addr);
	x = (double **) (op1->addr);
	y = (double **) (op2->addr);
	for (i = 0; i < (op1->dim)[0]; i++) {
		for (j = 0; j < (op1->dim)[1]; j++) {
			res[i][j] = x[i][j] + y[i][j];
		}
	}
}

void	vscaleinit(item *operator)
{
	item	*result, *op1, *op2;
	int	i;
	double	*res, *x, y;

	op2 = pop();
	op1 = pop();
	if ((op1->itype != NUM_TYPE) || (op2->itype != NUM_TYPE) ||
		(op1->ndim != 1) || (op2->ndim != 1) ||
		((op2->dim)[0] != 1)) {
		error_exit("vscale: bad arguments\n");
	}
	result = operator->next = newitem("VSCALE_TEMP", NUM_TYPE, 1, (op1->dim)[0]);
	res = (double *) calloc((op1->dim)[0], sizeof(double));
	result->addr = (void *) res;
	x = (double *) (op1->addr);
	y = *((double *) (op2->addr));
	for (i = 0; i < (op1->dim)[0]; i++) {
		res[i] = x[i] * y;
	}
}

void	vscaledoit(item *operator)
{
	item	*op1, *op2;
	int	i;
	double	*res, *x, y;

	op2 = pop();
	op1 = pop();
	res = (double *) ((operator->next)->addr);
	x = (double *) (op1->addr);
	y = *((double *) (op2->addr));
	for (i = 0; i < (op1->dim)[0]; i++) {
		res[i] = x[i] * y;
	}
}

void	mscaleinit(item *operator)
{
	item	*result, *op1, *op2;
	int	i, j;
	double	**res, **x, y;

	op2 = pop();
	op1 = pop();
	if ((op1->itype != NUM_TYPE) || (op2->itype != NUM_TYPE) ||
		(op1->ndim != 2) || (op2->ndim != 1) ||
		((op2->dim)[0] != 1)) {
		error_exit("mscale: bad arguments\n");
	}
	result = operator->next = newitem("MSCALE_TEMP", NUM_TYPE, 2, (op1->dim)[0], (op1->dim)[1]);
	res = (double **) calloc((op1->dim)[0], sizeof(double *));
	for (i = 0; i < (op1->dim)[0]; i++) {
		res[i] = (double *) calloc((op1->dim)[1], sizeof(double));
	}
	result->addr = (void *) res;
	x = (double * *) (op1->addr);
	y = *((double *) (op2->addr));
	for (i = 0; i < (op1->dim)[0]; i++) {
		for (j = 0; j < (op1->dim)[1]; j++) {
			res[i][j] = x[i][j] * y;
		}
	}
}

void	mscaledoit(item *operator)
{
	item	*op1, *op2;
	int	i, j;
	double	**res, **x, y;

	op2 = pop();
	op1 = pop();
	res = (double **) ((operator->next)->addr);
	x = (double **) (op1->addr);
	y = *((double *) (op2->addr));
	for (i = 0; i < (op1->dim)[0]; i++) {
		for (j = 0; j < (op1->dim)[1]; j++) {
			res[i][j] = x[i][j] * y;
		}
	}
}

void	vshiftinit(item *operator)
{
	item	*result, *op1, *op2;
	int	i;
	double	*res, *x, y;

	op2 = pop();
	op1 = pop();
	if ((op1->itype != NUM_TYPE) || (op2->itype != NUM_TYPE) ||
		(op1->ndim != 1) || (op2->ndim != 1) ||
		((op2->dim)[0] != 1)) {
		error_exit("vscale: bad arguments\n");
	}
	result = operator->next = newitem("VSCALE_TEMP", NUM_TYPE, 1, (op1->dim)[0]);
	res = (double *) calloc((op1->dim)[0], sizeof(double));
	result->addr = (void *) res;
	x = (double *) (op1->addr);
	y = *((double *) (op2->addr));
	for (i = 0; i < (op1->dim)[0]; i++) {
		res[i] = x[i] + y;
	}
}

void	vshiftdoit(item *operator)
{
	item	*op1, *op2;
	int	i;
	double	*res, *x, y;

	op2 = pop();
	op1 = pop();
	res = (double *) ((operator->next)->addr);
	x = (double *) (op1->addr);
	y = *((double *) (op2->addr));
	for (i = 0; i < (op1->dim)[0]; i++) {
		res[i] = x[i] + y;
	}
}


void	inverseinit(item *operator)
{
	item	*result, *op;
	int	i, j, dim;
	double	**res, **a, det;

	op = pop();
	if ((op->itype != NUM_TYPE) || (op->ndim != 2)) {
		error_exit("inverse: argument must be a 2-dimensional matrix\n");
	}
	dim = (op->dim)[0];
	if (dim != 2 && dim != 3) {
		error_exit("inverse: I can only handle matrices of size 2 or 3\n");
	}
	if ((op->dim)[1] != dim) {
		error_exit("inverse: matrix must be square\n");
	}
	result = operator->next = newitem("INVERSE_TEMP", NUM_TYPE, 2, dim, dim);
	res = (double **) calloc(dim, sizeof(double *));
	for (i = 0; i < dim; i++) {
		res[i] = (double *) calloc(dim, sizeof(double));
	}
	result->addr = (void *) res;
	a = (double **) (op->addr);
}

void	inversedoit(item *operator)
{
	item	*op;
	int	i, j, dim;
	double	**res, **a, det, invdet;

	op = pop();
	dim = (op->dim)[0];
	res = (double **) ((operator->next)->addr);
	a = (double **) (op->addr);
	if (dim == 2) {
		det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
		if (det != 0.0) {
			res[0][0] =  a[1][1] / det;
			res[0][1] = -a[0][1] / det;
			res[1][0] = -a[1][0] / det;
			res[1][1] =  a[0][0] / det;
		} else {
			res[0][0] = res[0][1] = res[1][0] = res[1][1] = 0.0;
		}
	} else {
		det = 	a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
			a[0][1] * (a[1][0] * a[2][2] - a[2][0] * a[1][2]) +
			a[0][2] * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
		invdet = (det == 0.0 ? 0.0 : 1.0 / det);
		res[0][0] =  invdet * (a[1][1] * a[2][2] - a[1][2] * a[2][1]);
		res[0][1] = -invdet * (a[0][1] * a[2][2] - a[2][1] * a[0][2]);	
		res[0][2] =  invdet * (a[0][1] * a[1][2] - a[1][1] * a[0][2]);
		res[1][0] = -invdet * (a[1][0] * a[2][2] - a[2][0] * a[1][2]);
		res[1][1] =  invdet * (a[0][0] * a[2][2] - a[2][0] * a[0][2]);
		res[1][2] = -invdet * (a[0][0] * a[1][2] - a[1][0] * a[0][2]);
		res[2][0] =  invdet * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
		res[2][1] = -invdet * (a[0][0] * a[2][1] - a[2][0] * a[0][1]);
		res[2][2] =  invdet * (a[0][0] * a[1][1] - a[1][0] * a[0][1]);	
	}
}