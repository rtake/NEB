# include <cstdio>
# include <cstdlib>
# include <cmath>


typedef double (*Func)(double arg);


typedef struct BaseInfo_{
	int argdim; // dimension of argument
	int* nexpvec; // vector of number of exponent
	Func* fvec; // vector of Function
}BaseInfo; // information of base function


double BaseFunc(BaseInfo* b, double *vec) {
	int i;
	double val = 1;
	for(int i = 0;i < b->argdim;i++) val *= pow(b->fvec[i](vec[i]),b->nexpvec[i]);
	return val;
} // BaseFunc include some Func and express as multiplication 


typedef struct Model_{
	int nbase; // number of base func
	int argdim; // dimension
	double* coeff; // coefficients
	BaseInfo** binfo; // vector of <information for base function>
	// Basefunc* bfvec; // vector of Base func
}Model;


double ModelFunc(Model *m, double* vec) {
	int i;
	double val = 0;
	for(i = 0;i < m->nbase;i++) { val += m->coeff[i] * BaseFunc(m->binfo[i],vec); }
	return val;
} // 


double d_ModelFunc(Model* m, double* xvec, double* dxvec) {
	int i;
	double df, dx;
	double* xxvec;

	dx = 0; for(int i = 0;i < m->argdim;i++) dx += dxvec[i] * dxvec[i]; dx = sqrt(dx); // set dx

	xxvec = (double*)malloc(sizeof(double) * m->argdim); // new xvec
	for(i = 0;i < m->argdim;i++) { xxvec[i] = xvec[i] + dxvec[i]; }
	df = ( ModelFunc(m,xxvec) - ModelFunc(m,xvec) ) / dx;
	free(xxvec);

	return df;
}


double Exptype(double r) { return 1 - exp(-0.5 * r);}


int Combination(int n, int k) { return 0; }


int** MakeCombination(int vecsize, int condition) {
	int i,j,k, matsize = Combination(vecsize + condition,condition);
	int* vec;
	int** imat; int** buf;
	imat = (int**)malloc(sizeof(int*) * matsize);
	for(i = 0;i < matsize;i++) imat[i] = (int*)malloc(sizeof(int) * vecsize);

	if(vecsize == 1) {
		for(i = 0;i < condition + 1;i++) {
			imat[i] = (int*)malloc(sizeof(int) * vecsize);
			for(j = 0;j < vecsize;j++) imat[i][j] = j;
		}
	} else {
		for(i = 0;i < condition + 1;i++) {
			buf = MakeCombination(vecsize - 1,i);


		}
	}

	return imat;
}


int main(int argc, char* argv[]) {
	int i,j,k;

	int argdim = 15, condition = 8, nbase = Combination(argdim + condition,condition);
	int** imat = MakeCombination(argdim,condition);

	Model *m;
	m = (Model*)malloc(sizeof(Model));
	m->nbase = nbase;
	m->argdim = argdim;
	m->coeff = (double*)malloc(sizeof(double) * nbase);
	m->binfo = (BaseInfo**)malloc(sizeof(BaseInfo*) * nbase);

	for(i = 0;i < nbase;i++) { // for each base function, ...
		// m->binfo->nexpvec = (int*)malloc(sizeof(int) * argdim);
		m->binfo[i]->nexpvec = imat[i];
		m->binfo[i]->fvec = (Func*)malloc(sizeof(Func) * argdim);
		for(j = 0;j < argdim;j++) { m->binfo[i]->fvec[j] = Exptype; }
	}




	for(i = 0;i < nbase;i++) free(m->binfo[i]->fvec);

	free(m->binfo);
	free(m->coeff);
	free(m);

	for(i = 0;i < nbase;i++) free(imat[i]);
	free(imat);

	return 0;
}
