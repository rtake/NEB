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

double BaseFunc(BaseInfo* b, double *vec) {
	int i;
	double val = 1;
	for(int i = 0;i < b->nexpvec;i++) val *= pow( b->fvec[i](vec[i]), b->expvec[i] );
	return val;
} //


double ModelFunc(Model* m, double* vec, int dim) {
	int i, j, k, npair, index;
	const int natom = dim / 3;
	double val, dist;
	double* argvec;
	
	npair = natom * (npair - 1);

	argvec = (double*)malloc(sizeof(double) * npair);

	index = 0;
	for(i = 0;i < natom;i++) {
		for(j = i + 1;j < natom;j++) {
			dist = 0; for(k = 0;k < 3;k++) { dist += (vec[i * 3 + k] - vec[j * 3 + k]) * (vec[i * 3 + k] - vec[j * 3 + k]); }
			argvec[index++] = sqrt(dist);
		}
	}

	for(i = 0;i < m->nbase;i++) { val += m->coeff[i] * BaseFunc(m->binfo[i],argvec); }

	free(argvec);

	return val;
} // Information for Model Potential // argument vec : 


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
		df[i] = ( ModelFunc(neb->m, vec_plus_dxvec, neb->argdim) - ModelFunc(neb->m, neb->images[num], neb->argdim) ) / (neb->d); // ok
	}
	
	free(vec_plus_dxvec);
	free(dxvec);

	return df;
} // Partial differential


double Exptype(double r) { return 1 - exp(-0.5 * r);}


void MakeIGuess(NEBInfo* neb) {
	int i, j, k;

	for(i = 1;i < neb->nimage - 1;i++) { // for each image except for the first and the last images
		for(j = 0;j < neb->argdim;j++) {
			neb->images[i][j] = i * (neb->images[neb->argdim - 1][j] - neb->images[0][j]) / neb->nimage;
		}
	}

}


double* Force_spring(NEBInfo* neb, int num) {
	int i;
	double* fvec;

	fvec = (double*)malloc(sizeof(double) * neb->argdim);

	for(i = 0;i < neb->argdim;i++) {
		fvec[i] = neb->k * ( (neb->images[num + 1][i] - neb->images[num][i]) - (neb->images[num][i] - neb->images[num - 1][i]) );
	}

	return fvec;
}


double* projection(int size, double* vec, double* t) {
	int i;
	double innerproduct;
	double* vec_projection;

	vec_projection = (double*)malloc(sizeof(double) * size);

	innerproduct = 0;
	for(i = 0;i < size;i++) { innerproduct += vec[i] * t[i]; }
	for(i = 0;i < size;i++) { vec_projection[i] = t[i] / innerproduct; }

	return vec_projection;
}


double** SetTanMat(NEBInfo* neb) {
	int i,j;

	double absvec;
	double* vec;
	double** tanmat; // set of tangent vector
	
	tanmat = (double**)malloc(sizeof(double*) * neb->nimage);
	for(i = 0;i < neb->nimage;i++) { tanmat[i] = (double*)malloc(sizeof(double) * neb->argdim); }

	for(i = 1;i < neb->nimage - 1;i++) { // for each images excepting for first and last images, ...

		absvec = 0;
		for(j = 0;j < neb->argdim;j++) {
			vec[j] = ( neb->images[i + 1] - neb->images[i - 1] );
			absvec += vec[j] * vec[j];
		}
		absvec = sqrt(absvec);

		for(j = 0;j < neb->argdim;j++) { tanmat[i][j] = vec[j] / absvec; }
	}

	free(vec);

	return tanmat;
} // return set of tangent vector


double* Force(NEBInfo* neb, int num) {
	int i;
	double* Fvec;
	double* Fs = Force_spring(neb,num); // malloc // ok
	double* Fs_projection = projection(neb->argdim, Fs, neb->tanmat[i]); // malloc // ok
	double* dV = d_ModelFunc(neb,num); // malloc // ok
	double* dV_projection = projection(neb->argdim, dV, neb->tanmat[i]); // malloc // ok

	Fvec = (double*)malloc(sizeof(double) * neb->argdim);
	for(i = 0;i < neb->argdim;i++) { Fvec[i] = Fs_projection[i] - dV_projection[i]; }

	free(dV); free(dV_projection);
	free(Fs); free(Fs_projection);
	return Fvec;
}

/*
double d_Force(NEBInfo* neb, double* dxvec) {
	double df;
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
		neb->tanmat = SetTanMat(neb); // malloc
		for(i = 1;i < neb->nimage - 1;i++) { // for each images
			dvec = Force(neb,i); // get gradient // malloc // ok 
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
	const int natom = 3, nimage = 10;
	const int argdim = natom * 3;

	int condition = 8, nbase = Combination(argdim + condition,condition), nimage = 10;
	int** imat = MakeCombination(argdim,condition);
	NEBInfo* nebinfo;



	nebinfo = (NEBInfo*)malloc(sizeof(NEBInfo));
	nebinfo->images = (double**)malloc(sizeof(double*) * nimage);
	for(i = 0;i < nimage;i++) { nebinfo->images[i] = (double*)malloc(sizeof(double) * argdim); }
	nebinfo->tanmat = (double**)malloc(sizeof(double*) * nimage);
	for(i = 0;i < nimage;i++) { nebinfo->images[i] = (double*)malloc(sizeof(double) * argdim); }

	nebinfo->m = (Model*)malloc(sizeof(Model));
	nebinfo->m->nbase = nbase;
	nebinfo->m->coeff = (double*)malloc(sizeof(double) * nbase);
	nebinfo->m->binfo = (BaseInfo**)malloc(sizeof(BaseInfo*) * nbase);

	for(i = 0;i < nbase;i++) { // for each base function, ...
		nebinfo->m->binfo[i]->expvec = imat[i];
		nebinfo->m->binfo[i]->fvec = (Func*)malloc(sizeof(Func) * nebinfo->m->nexpvec);
		for(j = 0;j < nexpvec;j++) { nebinfo->m->binfo[i]->fvec[j] = Exptype; }
	}

	// Memory Allocation END
	

	// Start NEB process

	MakeIGuess(nebinfo); // here

	Optimization_ModelFunc(nebinfo); // ok

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
