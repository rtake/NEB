# include <cstdio>
# include <cstdlib>
# include <cmath>
# include <cstring>


//////////////////////////////////////////////////////////////////////////////////////////////////

typedef double (*Func)(double arg);
// typedef double** Coordinate;


typedef struct BaseInfo_{
	int* expvec;
	int npair; // vector of number of exponent
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
	int natom;
	int order;
	int argdim;
	int npair;
	int nbase;

	char** elms;

	double** tanmat; // matrix of tangent
	int **imat;

	int k; // spring constant
	double d; // difference step
	double threshold;
	double alpha;

	char name[256];
	char log[256];
	char pes[256];
	char xyz[256]; // initial structure
}NEBInfo;

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

double BaseFunc(BaseInfo* b, double *vec) {
	int i;
	double val = 1;
	
	for(i = 0;i < b->npair;i++) val *= pow( b->fvec[i](vec[i]), b->expvec[i] ); // here
	// printf("val %lf",val); for(i = 0;i < b->npair;i++) printf(" %lf",vec[i]); printf("\n");
	
	return val;
} //


double ModelFunc(Model* m, double* vec, int dim) {
	int i, j, k, npair, index;
	const int natom = dim / 3;
	double val, dist;
	double* argvec;
	
	npair = natom * (natom - 1);
	argvec = (double*)malloc(sizeof(double) * npair);

	index = 0;
	for(i = 0;i < natom;i++) { // for each atom, ...
		for(j = i + 1;j < natom;j++) { // for each atom, ...
			dist = 0; for(k = 0;k < 3;k++) { dist += (vec[3*i + k] - vec[3*j + k]) * (vec[3*i + k] - vec[3*j + k]); }
			argvec[index] = sqrt(dist);
			// printf("dist %d - %d : %lf\n",i,j,argvec[index]);
			index++;
		}
	}

	val = 0;
	for(i = 0;i < m->nbase;i++) { val += m->coeff[i] * BaseFunc(m->binfo[i],argvec); } // here

	free(argvec);

	return val;
} // Information for Model Potential // argument vec : 


double* d_ModelFunc(NEBInfo* neb, int num) { // for neb->images[num]
	int i,j;
	double* df; double* dxvec; double* vec_plus_dxvec;

	df = (double*)malloc(sizeof(double) * neb->argdim);
	dxvec = (double*)malloc(sizeof(double) * neb->argdim);
	vec_plus_dxvec = (double*)malloc(sizeof(double) * neb->argdim);

	for(i = 0;i < neb->argdim;i++) { // for each df elems, ...
		
		for(j = 0;j < neb->argdim;j++) {
			if(j == i) dxvec[j] = neb->d;
			else dxvec[j] = 0;
		}
		
		for(j = 0;j < neb->argdim;j++) { vec_plus_dxvec[j] = neb->images[num][j] + dxvec[j]; }
		
		df[i] = ( ModelFunc(neb->m, vec_plus_dxvec, neb->argdim) - ModelFunc(neb->m, neb->images[num], neb->argdim) ) / (neb->d); // here
		
		printf("vec_plus_dxvec   :"); for(j = 0;j < neb->argdim;j++) printf(" %lf",vec_plus_dxvec[j]); printf("\n");
		printf("neb->images[%d]  :",num); for(j = 0;j < neb->argdim;j++) printf(" %lf",neb->images[num][j]); printf("\n");
	}

	
	printf("df:\n"); for(i = 0;i < neb->argdim;i++) printf("%lf\t",df[i]); printf("\n");

	free(vec_plus_dxvec);
	free(dxvec);

	return df;
} // Partial differential


double Exptype(double r) { return 1 - exp(-0.5 * r);}


void LoadIniStr(NEBInfo* neb, FILE* fp) {
	int i, j, natom;
	char buf[256], line[256];
	
	for(i = 0;i < 2;i++) {
		fgets(line,256,fp);
		sscanf(line,"%d",&natom);
		fgets(line,256,fp);
		for(j = 0;j < natom;j++) {
			fgets(line,256,fp); // printf("%s",line);
			if(i == 0) sscanf(line,"%s%17lf%17lf%17lf",neb->elms[j],&neb->images[0][3*j],&neb->images[0][3*j + 1],&neb->images[0][3*j + 2]);
			else if(i == 1) sscanf(line,"%s%17lf%17lf%17lf",buf,&neb->images[neb->nimage - 1][3*j],&neb->images[neb->nimage - 1][3*j + 1],&neb->images[neb->nimage - 1][3*j + 2]);
		}
	}

} // load initial structures ( first, end )


void MakeIGuess(NEBInfo* neb) {
	int i, j, k;
	for(i = 1;i < neb->nimage - 1;i++) { // for each image except for the first and the last images
		for(j = 0;j < neb->argdim;j++) {
			neb->images[i][j] = neb->images[0][j] + i * (neb->images[neb->nimage - 1][j] - neb->images[0][j]) / (neb->nimage - 1);
			// printf("images[%d][%d] : %lf\t", i, j, neb->images[i][j]);
		} // printf("\n");
	}

	printf("Initial Guess\n");
	for(i = 0;i < neb->nimage;i++) {
		printf("%d\n%d\n",neb->argdim ,i);
		for(j = 0;j < neb->argdim / 3;j++) { printf("%s\t%17lf\t%17lf\t%17lf\n",neb->elms[j], neb->images[i][3*j], neb->images[i][3*j + 1], neb->images[i][3*j + 2]); }
	}	

}


double* Force_spring(NEBInfo* neb, int num) {
	int i;
	double* fvec;

	fvec = (double*)malloc(sizeof(double) * neb->argdim);

	for(i = 0;i < neb->argdim;i++) {
		fvec[i] = neb->k * ( (neb->images[num + 1][i] - neb->images[num][i]) - (neb->images[num][i] - neb->images[num - 1][i]) );
		// printf("images ([%d][%d] - [%d][%d]) : %lf, ([%d][%d] - [%d][%d]) : %lf\n",num + 1, i, num, i, (neb->images[num + 1][i] - neb->images[num][i]), num, i, num - 1, i, (neb->images[num][i] - neb->images[num - 1][i]));
	}

	return fvec;
}


double* projection(int size, double* vec, double* t) {
	int i;
	double innerproduct;
	double* vec_projection;

	vec_projection = (double*)malloc(sizeof(double) * size);

	innerproduct = 0;

	for(i = 0;i < size;i++) printf("vec[%d] %lf, t[%d] %lf\t",i, vec[i], i, t[i]); printf("\n");
	for(i = 0;i < size;i++) innerproduct += vec[i] * t[i]; printf("innerproduct %lf\n",innerproduct);
	
	for(i = 0;i < size;i++) vec_projection[i] = t[i] * innerproduct;
	// printf("vec_projection:\n"); for(i = 0;i < size;i++) printf("%lf\t",vec_projection[i]);
	
	return vec_projection;
}


void SetTanMat(NEBInfo* neb) {
	int i,j;
	double absvec;
	double* vec;

	vec = (double*)malloc(sizeof(double) * neb->argdim);
	for(i = 1;i < neb->nimage - 1;i++) { // for each images excepting for first and last images, ...
		absvec = 0;
		for(j = 0;j < neb->argdim;j++) {
			vec[j] = ( neb->images[i + 1][j] - neb->images[i - 1][j] );
			absvec += vec[j] * vec[j];
			// printf("%lf = ( neb->images[%d + 1][%d] - neb->images[%d - 1][%d] )\n",vec[j], i, j, i, j);
		}
		absvec = sqrt(absvec); // printf("absvec %lf\n",absvec);
		for(j = 0;j < neb->argdim;j++) { neb->tanmat[i][j] = vec[j] / absvec; }
	}
	free(vec);

	printf("tanmat\n");
	for(i = 1;i < neb->nimage - 1;i++) { 
		for(j = 0;j < neb->argdim;j++) {
			printf("%d %d %lf\t",i,j,neb->tanmat[i][j]);
		}
		printf("\n");
	}

} // return set of tangent vector


double* Force(NEBInfo* neb, int num) { // printf("Force %d start\n",num);
	int i;
	double* Fvec;
	double* Fs = Force_spring(neb,num); // malloc // ok
	double* Fs_projection = projection(neb->argdim, Fs, neb->tanmat[num]); // malloc // ok
	double* dV = d_ModelFunc(neb,num); // malloc // here
	double* dV_projection = projection(neb->argdim, dV, neb->tanmat[num]); // malloc // ok

	Fvec = (double*)malloc(sizeof(double) * neb->argdim);
	for(i = 0;i < neb->argdim;i++) Fvec[i] = Fs_projection[i] - dV_projection[i]; // hehe NAN 
	for(i = 0;i < neb->argdim;i++) {
		printf("Fs[%d] (%lf), dV[%d] (%lf)\n",i, Fs[i], i, dV[i]);
		printf("Fvec[%d] (%lf)  = Fs_projection[%d] (%lf) - dV_projection[%d] (%lf)\n",i, Fvec[i], i, Fs_projection[i], i, dV_projection[i]);
	}

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
	int i, j, chk = 0, cycle;
	const int maxcycle = 100;
	double *dvec = NULL, **dmat = NULL; // grad vector
	FILE *fp;
	char fname[256], name[256];

	dmat = (double**)malloc(sizeof(double*) * neb->nimage);
	for(i = 0;i < neb->nimage;i++) dmat[i] = (double*)malloc(sizeof(double) * neb->argdim);

	for(cycle = 0;cycle < maxcycle;cycle++) { printf("OptCycle %d start\n",cycle);
		chk = 1; // check converged
		SetTanMat(neb); // ok

		for(i = 1;i < neb->nimage - 1;i++) { // for each images
			dvec = Force(neb,i); // get gradient // malloc
			for(j = 0;j < neb->argdim;j++) dmat[i][j] = dvec[j];
			for(j = 0;j < neb->argdim;j++) { if(dvec[j] > neb->threshold) chk = 0; }
			free(dvec);
		}

		printf("dmat(cycle %d):\n",cycle);
		for(i = 0;i < neb->nimage;i++) {
			for(j = 0;j < neb->argdim;j++) {
				if(i == 0 || i == neb->nimage - 1) printf("0 ");
				else printf("%lf ",dmat[i][j]);
			} printf("\n");
		}

		if(chk == 1) break; // converged

		for(i = 1;i < neb->nimage - 1;i++) { 
			for(j = 0;j < neb->argdim;j++) { neb->images[i][j] -= (neb->alpha) * dmat[i][j]; }
		} // update coordinate

		if(cycle %10 == 0) {
			sprintf(fname,"%s_NEB_cycle%d.neb.xyz", neb->name, cycle);
			fp = fopen(fname,"w");
			for(i = 0;i < neb->nimage;i++) {
				fprintf(fp,"%d\n%d\n",neb->argdim / 3, i);
				for(j = 0;j < neb->argdim / 3;j++) { // printf("i %d j %d\n",i,j);
					fprintf(fp,"%s\t%17lf\t%17lf\t%17lf\n",neb->elms[j], neb->images[i][3*j], neb->images[i][3*j + 1], neb->images[i][3*j + 2]);
				}
			}
			fclose(fp); // FILE OUT
		} // FILE OUT

	} // Optimization (steepest decent)

	for(i = 0;i < neb->argdim;i++) free(dmat[i]); free(dmat);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

int Combination(int n, int r) {
	if(n == r) return 1;
	else if(r == 0) return 1;
	else if(r == 1) return n;
	else return Combination(n - 1, r- 1) + Combination(n - 1, r);
}


int** MakeCombination(int vecsize, int condition) { // printf("vecsize %d, condition %d\n",vecsize, condition);
	int index, i,j,k, matsize_buf, matsize = Combination(vecsize + condition,condition); // printf("matsize_buf : %d\n",matsize);
	int* vec;
	int** imat; int** buf;
	imat = (int**)malloc(sizeof(int*) * matsize);
	for(i = 0;i < matsize;i++) imat[i] = (int*)malloc(sizeof(int) * vecsize);

	index = 0;
	if(vecsize == 1) {
		for(i = 0;i < condition + 1;i++) {
			imat[i][index] = i; // printf("vecsize 1 imat[%d][0] : %d\n",i, imat[i][0]);
		}
	} else {
		for(i = 0;i < condition + 1;i++) { // for each condition
			buf = MakeCombination(vecsize - 1,i);
			matsize_buf = Combination(vecsize - 1 + i,i);

			for(j = 0;j < matsize_buf;j++, index++) { // for each row, ...
				imat[index][0] = condition - i;
				for(k = 0;k < vecsize - 1;k++) {
					// printf("buf[%d][%d] : %d\t",j, k, buf[j][k]);
					imat[index][k + 1] = buf[j][k];
				}
				// printf("\n");
			}
		}
	}

/*
	printf("imat\n");
	for(i = 0;i < matsize;i++) { for(j = 0;j < vecsize;j++) { printf("%d ",imat[i][j]); } printf("\n"); }
	printf("vecsize %d, condition %d end\n\n",vecsize, condition);
*/

	return imat;
}


void LoadInfile(NEBInfo* neb) {
	FILE *fp;
	char line[256], *pt, infname[256];

	fp = fopen( infname, "r" );
	while( fgets( line, 256, fp ) ) {
		pt = strstr( line, ":" );
		if( strstr( line, "pes") ) sscanf( pt + 2, "%s", neb->pes );
		else if( strstr( line, "ini") ) sscanf( pt + 2, "%s", neb->xyz );
		else if( strstr( line, "log") ) sscanf( pt + 2, "%s", neb->log );
		else if( strstr( line, "image") ) sscanf( pt + 2, "%d", &neb->nimage );
		else if( strstr( line, "atom") ) sscanf( pt + 2, "%s", &neb->natom );
		else if( strstr( line, "order") ) sscanf( pt + 2, "%s", &neb->order );
	}
	fclose( fp );

	neb->argdim = neb->natom * 3;
	neb->npair = neb->natom * ( neb->natom - 1 ) / 2;
	neb->nbase = Combination( neb->npair + neb->order, neb->order );
}

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
	int i,j,k;
	char infname[256], fname[256], xyz[256], neb[256], line[256], pes[256], *pt;
	FILE *fp_xyz, *fp_neb, *fp_pes, *fp_in;
	NEBInfo* nebinfo;

	nebinfo = (NEBInfo*)malloc( sizeof(NEBInfo) ); // resore calc. info.

	sprintf( nebinfo->name, "%s", argv[1] );
	LoadInfile( nebinfo );

	nebinfo->imat = MakeCombination( nebinfo->npair, nebinfo->order );
	nebinfo->k = 1.0;
	nebinfo->d = 0.0001;
	nebinfo->threshold = 0.0001;
	nebinfo->alpha = 0.1;

	// Memory Allocation

	nebinfo->images = (double**)malloc( sizeof(double*) * nebinfo->nimage );
	for(i = 0;i < nebinfo->nimage;i++) { nebinfo->images[i] = (double*)malloc(sizeof(double) * nebinfo->argdim); }
	nebinfo->tanmat = (double**)malloc(sizeof(double*) * nebinfo->nimage);
	for(i = 0;i < nebinfo->nimage;i++) { nebinfo->tanmat[i] = (double*)malloc(sizeof(double) * nebinfo->argdim); }

	nebinfo->elms = (char**)malloc(sizeof(char*) * nebinfo->natom);
	for(i = 0;i < nebinfo->natom;i++) { nebinfo->elms[i] = (char*)malloc(sizeof(char) * 10); }

	nebinfo->m = (Model*)malloc(sizeof(Model));
	nebinfo->m->coeff = (double*)malloc(sizeof(double) * nebinfo->nbase);
	nebinfo->m->binfo = (BaseInfo**)malloc(sizeof(BaseInfo*) * nebinfo->nbase);

	for(i = 0;i < nebinfo->nbase;i++) { // for each base function, ...
		nebinfo->m->binfo[i] = (BaseInfo*)malloc(sizeof(BaseInfo));
		nebinfo->m->binfo[i]->expvec = nebinfo->imat[i];
		nebinfo->m->binfo[i]->npair = nebinfo->npair;
		nebinfo->m->binfo[i]->fvec = (Func*)malloc(sizeof(Func) * nebinfo->npair);
		for(j = 0;j < nebinfo->npair;j++) { nebinfo->m->binfo[i]->fvec[j] = Exptype; }
	}

	// Memory Allocation END

	// Start NEB process

	LoadIniStr( nebinfo );

	LoadPES( nebinfo );

	/*
	fp_pes = fopen( pes, "r" );
	for(i = 0;i < nebinfo->nbase;i++) {
		fgets( line, 256, fp_pes );
		sscanf( line, "%lf", &nebinfo->m->coeff[i] );
	}
	fclose(fp_pes);
	*/

	MakeIGuess( nebinfo ); // ok

	Optimization_ModelFunc( nebinfo ); // here

	// End NEB process

	// FILE Out 

	fp_neb = fopen(neb,"w");
	for(i = 0;i < nebinfo->nimage;i++) {
		fprintf(fp_neb,"%d\n%d\n", nebinfo->natom, i);
		for(j = 0;j < nebinfo->natom;j++) {
			fprintf(fp_neb,"%s\t%17lf\t%17lf\t%17lf\n", nebinfo->elms[j], nebinfo->images[i][3*j], nebinfo->images[i][3*j + 1], nebinfo->images[i][3*j + 2]);
		}
	}
	fclose(fp_neb);

	// FILE Out END

	// Free Memory

	for(i = 0;i < nebinfo->nbase;i++) free(nebinfo->m->binfo[i]->fvec);
	free(nebinfo->m->binfo);
	free(nebinfo->m->coeff);
	free(nebinfo->m);

	for(i = 0;i < nebinfo->natom;i++) { free(nebinfo->elms[i]); } free(nebinfo->elms);
	for(i = 0;i < nebinfo->nimage;i++) { free(nebinfo->tanmat[i]); } free(nebinfo->tanmat);
	for(i = 0;i < nebinfo->nimage;i++) { free(nebinfo->images[i]); } free(nebinfo->images);

	free( nebinfo->imat );
	free(nebinfo);

	// Free Memory END

	return 0;
}
