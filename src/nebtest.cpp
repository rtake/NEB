# include <cstdio>
# include <cstdlib>
# include <cmath>


//////////////////////////////////////////////////////////////////////////////////////////////////

typedef double (*Func)(double arg);
// typedef double** Coordinate;


typedef struct BaseInfo_{
	int* expvec;
	int nexpvec; // vector of number of exponent
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

	double** images; // coordinate vector
	int nimage; // number of images
	int argdim;

	double** tanmat; // matrix of tangent

	int k; // spring constant

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
	for(int i = 0;i < b->nexpvec;i++) val *= pow( b->fvec[i](vec[i]), b->expvec[i] );
	return val;
} //


double ModelFunc(Model *m, double* vec) {
	int i;
	double val = 0;
	// for(i = 0;i < m->nbase;i++) { val += m->coeff[i] * BaseFunc(m->binfo[i],vec); } // here
	return val;
} // Information for Model Potential


double* d_ModelFunc(NEBInfo* neb, int num) {
	int i,j;
	double* df; double* dxvec; double* vec_plus_dxvec;

	dxvec = (double*)malloc(sizeof(double) * neb->argdim);
	vec_plus_dxvec = (double*)malloc(sizeof(double) * neb->argdim);

	for(i = 0;i < neb->argdim;i++) { // for each df elems, ...
		
		for(j = 0;j < neb->argdim;j++) {
			if(j == i) dxvec[j] = neb->d;
			else dxvec[j] = 0;
		}

		for(j = 0;j < neb->argdim;j++) { vec_plus_dxvec[j] = neb->images[num][j] + dxvec[j]; }
		df[i] = ( ModelFunc(neb->m,vec_plus_dxvec) - ModelFunc(neb->m,neb->images[num]) ) / (neb->d); // here
	}
	
	free(vec_plus_dxvec);
	free(dxvec);

	return df;
} // Partial differential


double Exptype(double r) { return 1 - exp(-0.5 * r);}


void MakeIGuess(NEBInfo* neb) {}


double* Force_spring(NEBInfo* neb, int num) {
	int i;
	double* fvec;

	return fvec;
}


double* projection(double* vec, double* t) {
	double* vec_projection;


	return vec_projection;
}


double** SetTanMat(NEBInfo* neb) {
	double** tanmat;


	return tanmat;
}


double* Force(NEBInfo* neb, int num) {
	int i;
	double* Fvec;
	double* Fs = Force_spring(neb,num); // malloc here
	double* Fs_projection = projection(Fs,neb->tanmat[i]); // malloc here
	double* dV = d_ModelFunc(neb,num); // malloc here
	double* dV_projection = projection(dV,neb->tanmat[i]); // malloc here

	Fvec = (double*)malloc(sizeof(double) * neb->argdim);
	for(i = 0;i < neb->argdim;i++) { Fvec[i] = Fs_projection[i] - dV_projection[i]; }

	free(dV); free(dV_projection);
	free(Fs); free(Fs_projection);
	return Fvec;
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
	double* dvec = NULL; // grad vector
	// double** dxmat; // vector of (dx vector)
	
	// dvec = (double*)malloc(sizeof(double) * neb->argdim); // vector of difference
	// dxmat = (double**)malloc(sizeof(double) * neb->argdim);

	/*
	for(i = 0;i < neb->argdim;i++) {
		dxmat[i] = (double*)malloc(sizeof(double) * neb->argdim);
		for(j = 0;j < neb->argdim;j++) {
			if(j == i) dxmat[i][j] = neb->d;
			else dxmat[i][j] = 0;
		}
	}
	*/

	while(chk == 0) {
		chk = 1; // check converged
		neb->tanmat = SetTanMat(neb); // malloc here
		for(i = 1;i < neb->nimage - 1;i++) { // for each images
			// dvec = Force(neb,i); // get gradient / malloc here
			for(j = 0;j < neb->argdim;j++) { if(dvec[j] > neb->threshold) chk = 0; }
			if(chk == 1) break; // converged
			for(j = 0;j < neb->argdim;j++) { neb->images[i][j] -= (neb->alpha) * dvec[j]; }
			free(dvec);
		} // update coordinate
		free(neb->tanmat);
	} // Optimization (steepest decent)

	// for(i = 0;i < neb->argdim;i++) free(dxmat[i]);
	// free(dxmat);
}

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

	/*
	for(i = 0;i < nebinfo->m->nbase;i++) { // for each base function, ...
		nebinfo->m->binfo[i]->expvec = imat[i];
		nebinfo->m->binfo[i]->fvec = (Func*)malloc(sizeof(Func) * nebinfo->m->nexpvec);
		for(j = 0;j < nexpvec;j++) { nebinfo->m->binfo[i]->fvec[j] = Exptype; }
	}
	*/

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
