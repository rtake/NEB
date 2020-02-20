# include <cstdio>
# include <cstdlib>
# include <cmath>


//////////////////////////////////////////////////////////////////////////////////////////////////

typedef double (*Func)(double arg);
// typedef double** Coordinate;


typedef struct BaseInfo_{
	// int argdim; // dimension of argument
	int* nexpvec; // vector of number of exponent
	Func* fvec; // vector of Function
}BaseInfo; // information of base function


typedef struct Model_{
	int nbase; // number of base func
	// int argdim; // dimension
	double* coeff; // coefficients
	BaseInfo** binfo; // vector of <information for base function>
}Model;


typedef struct NEBInfo_{
	Model* m; // Information for Model Function

	double* vec; // coordinate vector
	int argdim;

	int p; // number of Images
	double d; // difference step
	double threshold;
	double alpha;
}NEBInfo;

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

/*
double* InnerProduct(double* vec0, double* vec1, int n) {
	int i;
	double val = 0;
	for(i = 0;i < n;i++) val += vec0[i] * vec1[i];
	return val;
}
*/

double BaseFunc(BaseInfo* b, double *vec) {
	int i;
	double val = 1;
	for(int i = 0;i < b->argdim;i++) val *= pow(b->fvec[i](vec[i]),b->nexpvec[i]);
	return val;
} //


double ModelFunc(Model *m, double* vec) {
	int i;
	double val = 0;
	for(i = 0;i < m->nbase;i++) { val += m->coeff[i] * BaseFunc(m->binfo[i],vec); }
	return val;
} // Information for Model Potential


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
} // Partial differential


double Exptype(double r) { return 1 - exp(-0.5 * r);}


void MakeIGuess(NEBInfo* neb) {}


double* projection(double* vec, double* tvec) {
	double* vec_projection;


	return vec_projection;
}


// void Force(NEBInfo* neb, double* fvec) {
double* Force(NEBInfo* neb) {
	int i;
	double** unit_tan_mat;
	double* fvec;

	fvec = (double*)malloc(sizeof(double) * neb->argdim);

	for(i = 0;i < neb->argdim;i++) {
		fvec[i] = 0;
		fvec[i] += projection( Force_spring(), unit_tan_mat[i] ); // spring
		fvec[i] -= d_ModelFunc() - projection( d_ModelFunc(), unit_tan_mat[i] );
	}

	fvec;
}

/*
double d_Force(NEBInfo* neb, double* dxvec) {
	int i;
	double df, dx;
	double* vec_plus_dxvec;
	vec_plus_dxvec = (double*)malloc(sizeof(double) * neb->argdim); // new xvec

	dx = 0; for(int i = 0;i < neb->argdim;i++) { dx += dxvec[i] * dxvec[i]; }  dx = sqrt(dx);
	for(i = 0;i < neb->argdim;i++) { vec_plus_dxvec[i] = neb->vec[i] + dxvec[i]; }
	// df = ( Force() - Force() ) / dx; // here

	free(vec_plus_dxvec);
	return df; 
}
*/

void Optimization_ModelFunc(NEBInfo* neb) {
	int i, j, chk = 0;
	double* dvec; // grad vector
	double** dxmat; // vector of (dx vector)
	dvec = (double*)malloc(sizeof(double) * neb->argdim); // vector of difference
	dxmat = (double**)malloc(sizeof(double) * neb->argdim);
	for(i = 0;i < neb->argdim;i++) {
		dxmat[i] = (double*)malloc(sizeof(double) * neb->argdim);
		for(j = 0;j < neb->argdim;j++) {
			if(j == i) dxmat[i][j] = neb->d;
			else dxmat[i][j] = 0;
		}
	}

	while(chk == 0) {
		chk = 1;

		// Force(neb,dvec);
		dvec = Force(neb);

		/*
		for(int i = 0;i < neb->argdim;i++) {
			dvec[i] = d_Force(neb,dxmat[i]); // get gradient // here
			if(dvec[i] > neb->threshold) chk = 0; // continue while all elems converged
		} // differentiate & chk
		*/

		if(chk == 1) break; // converged

		for(int i = 0;i < neb->argdim;i++) { neb->vec[i] -= (neb->alpha) * dvec[i]; } // update
	} // Optimization

	for(i = 0;i < neb->argdim;i++) free(dxmat[i]);
	free(dxmat);
	free(dvec);
} // steepest decent

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

int Combination(int n, int k) { return 0; }


int** MakeCombination(int vecsize, int condition) {
	int i,j,k, matsize = Combination(vecsize + condition,condition);
	int* vec;
	int** imat; int** buf;
	imat = (int**)malloc(sizeof(int*) * matsize);
	for(i = 0;i < matsize;i++) imat[i] = (int*)malloc(sizeof(int) * vecsize);

	if(vecsize == 1) {

		for(i = 0;i < condition + 1;i++) {
			for(j = 0;j < vecsize;j++) { imat[i][j] = j; }
		}

	} else {

		for(i = 0;i < condition + 1;i++) {
			buf = MakeCombination(vecsize - 1,i);
			
			for(j = 0;j < matsize;j++) { // for each row, ...
				
				imat[j][0] = condition - i;
				for(k = 0;k < vecsize - 1;k++) { imat[j][k + 1] = buf[j][k]; }

			}

		}
	}

	return imat;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
	int i,j,k;

	// Memory Allocation

	NEBInfo* nebinfo;
	nebinfo = (NEBInfo*)malloc(sizeof(NEBInfo));

	int argdim = 15, condition = 8, nbase = Combination(argdim + condition,condition), nimage = 10;
	int** imat = MakeCombination(argdim,condition);

	nebinfo->m = (Model*)malloc(sizeof(Model));
	nebinfo->m->nbase = nbase;
	// nebinfo->m->argdim = argdim;
	nebinfo->m->coeff = (double*)malloc(sizeof(double) * nbase);
	nebinfo->m->binfo = (BaseInfo**)malloc(sizeof(BaseInfo*) * nbase);

	for(i = 0;i < nbase;i++) { // for each base function, ...
		nebinfo->m->binfo[i]->nexpvec = imat[i];
		nebinfo->m->binfo[i]->fvec = (Func*)malloc(sizeof(Func) * argdim);
		for(j = 0;j < argdim;j++) { nebinfo->m->binfo[i]->fvec[j] = Exptype; }
	}

	// Memory Allocation END
	

	// Start NEB process

	MakeIGuess(nebinfo); // here

	Optimization_ModelFunc(nebinfo); // here

	// End NEB process


	// Free Memory

	for(i = 0;i < nbase;i++) free(nebinfo->m->binfo[i]->fvec);

	free(nebinfo->m->binfo);
	free(nebinfo->m->coeff);
	free(nebinfo->m);

	for(i = 0;i < nbase;i++) free(imat[i]);
	free(imat);

	free(nebinfo);

	// Free Memory END

	return 0;
}
