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


typedef struct {
	Model* m; // Information for Model Function

	double** images; // coordinate vector
	int nimage; // number of images
	int natom;
	int order;
	int argdim;
	int npair;
	int nbase;
	int cycle;

	char** elms;

	double** tanmat; // matrix of tangent
	int **imat;

	double k; // spring constant
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
	
	for(i = 0;i < b->npair;i++) { val *= pow( b->fvec[i](vec[i]), b->expvec[i] ); }

	/* chk */
	/*
	fprintf( stdout, "val %lf", val );
	for(i = 0;i < b->npair;i++) { fprintf( stdout, " %lf", vec[i] ); }
	fprintf( stdout, "\n" );
	*/

	return val;
}


double ModelFunc(Model* m, double* vec, int dim) {
	int i, j, k, npair, index;
	const int natom = dim / 3;
	double val, dist;
	double* argvec;
	
	npair = natom * (natom - 1);
	argvec = (double*)malloc(sizeof(double) * npair);

	index = 0;
	for(i = 0;i < natom;i++) { // for each atom, ...
		for(j = i + 1;j < natom;j++, index++) { // for each atom, ...
			dist = 0; for(k = 0;k < 3;k++) { dist += (vec[3*i + k] - vec[3*j + k]) * (vec[3*i + k] - vec[3*j + k]); }
			argvec[index] = sqrt(dist);
			// fprintf( stdout, "dist %d - %d : %lf\n", i, j, argvec[index] );
		}
	}

	val = 0;
	for(i = 0;i < m->nbase;i++) { val += m->coeff[i] * BaseFunc(m->binfo[i],argvec); } // here

	free(argvec);

	return val;
} // Information for Model Potential // argument vec : 


double* d_ModelFunc(NEBInfo* neb, int num) { // for neb->images[num]
	int i,j;
	double *df, *dxvec, *vec_plus_dxvec, f_vec_plus_dxvec, f_vec;

	df = (double*)malloc(sizeof(double) * neb->argdim);
	dxvec = (double*)malloc(sizeof(double) * neb->argdim);
	vec_plus_dxvec = (double*)malloc(sizeof(double) * neb->argdim);

	for(i = 0;i < neb->argdim;i++) { // for each df elems, ...
		
		for(j = 0;j < neb->argdim;j++) {
			if(j == i) dxvec[j] = neb->d;
			else dxvec[j] = 0;
		} // set dxvec
		
		for(j = 0;j < neb->argdim;j++) { vec_plus_dxvec[j] = neb->images[num][j] + dxvec[j]; }
		
		f_vec_plus_dxvec = ModelFunc( neb->m, vec_plus_dxvec, neb->argdim ); // here
		f_vec = ModelFunc( neb->m, neb->images[num], neb->argdim ); // here
		df[i] = ( f_vec_plus_dxvec - f_vec ) / ( neb->d );
	
		/* chk */
		/*
		fprintf( stdout, "df[%d] : ( %lf - %lf ) / %lf = %lf\n", i, f_vec_plus_dxvec, f_vec, neb->d, df[i] );
		*/
		/* chk ok */
	}

	free(vec_plus_dxvec);
	free(dxvec);

	return df;
} // Partial differential


double Exptype(double r) { return 1 - exp(-0.5 * r);}


void LoadIniStr( NEBInfo* neb ) {
	int i, j;
	char buf[256], line[256];
	FILE *fp;

	fp = fopen( neb->xyz, "r" );
	for(i = 0;i < 2;i++) {
		fgets(line,256,fp); fgets(line,256,fp);
		for(j = 0;j < neb->natom;j++) {
			fgets(line,256,fp); // printf("%s",line);
			if(i == 0) sscanf(line,"%s%17lf%17lf%17lf",neb->elms[j],&neb->images[0][3*j],&neb->images[0][3*j + 1],&neb->images[0][3*j + 2]);
			else if(i == 1) sscanf( line, "%s%17lf%17lf%17lf", buf, &neb->images[neb->nimage - 1][3*j],&neb->images[neb->nimage - 1][3*j + 1],&neb->images[neb->nimage - 1][3*j + 2]);
		}
	}
	fclose(fp);

	/* chk */
	/*	
	for(i = 0;i < neb->nimage;i++) {
		printf( "%d\n%d\n", neb->natom, i );
		for(j = 0;j < neb->natom;j++) {
			printf( "%s%17.12lf%17.12lf%17.12lf\n", neb->elms[j], neb->images[i][3*j], neb->images[i][3*j + 1], neb->images[i][3*j + 2] );
		}
	}
	*/

} // load initial structures ( first, end )


void MakeIGuess(NEBInfo* neb) {
	int i, j;
	for(i = 1;i < neb->nimage - 1;i++) { // for each image except for the first and the last images
		for(j = 0;j < neb->argdim;j++) {
			neb->images[i][j] = neb->images[0][j] + i * (neb->images[neb->nimage - 1][j] - neb->images[0][j]) / (neb->nimage - 1);
			// printf("images[%d][%d] : %lf\t", i, j, neb->images[i][j]);
		} // printf("\n");
	}

	/* chk */
	/*
	printf("Initial Guess\n");
	for(i = 0;i < neb->nimage;i++) {
		fprintf( stdout, "%d\n%d\n", neb->natom, i );
		for(j = 0;j < neb->argdim / 3;j++) { fprintf( stdout, "%s\t%17.12lf\t%17.12lf\t%17.12lf\n",neb->elms[j], neb->images[i][3*j], neb->images[i][3*j + 1], neb->images[i][3*j + 2]); }
	}	
	*/
}


double* Force_spring(NEBInfo* neb, int num) {
	int i;
	double* fvec;

	fvec = (double*)malloc(sizeof(double) * neb->argdim);
	for(i = 0;i < neb->argdim;i++) {
		fvec[i] = neb->k * ( (neb->images[num + 1][i] - neb->images[num][i]) - (neb->images[num][i] - neb->images[num - 1][i]) );

		/* chk */
		/*
		fprintf( stdout, "Fs[%d][%d] : %lf\t", num, i, fvec[i] );
		fprintf( stdout, "([%d][%d] - [%d][%d]) : %lf\t", num + 1, i, num, i, ( neb->images[num + 1][i] - neb->images[num][i] ) );
		fprintf( stdout, "([%d][%d] - [%d][%d]) : %lf\n", num, i, num - 1, i, ( neb->images[num][i] - neb->images[num - 1][i] ) );
		*/
		/* ok */
	}

	return fvec;
}


double* projection(int size, double* vec, double* t) {
	int i;
	double innerproduct = 0;
	double* vec_projection;

	vec_projection = (double*)malloc(sizeof(double) * size);

	for(i = 0;i < size;i++) { innerproduct += vec[i] * t[i]; }
	for(i = 0;i < size;i++) vec_projection[i] = t[i] * innerproduct;

	/* chk */
	/*
	for(i = 0;i < size;i++) { fprintf( stdout, "vec[%d] %lf, t[%d] %lf\n", i, vec[i], i, t[i] ); }
	fprintf( stdout, "innerproduct %lf\n", innerproduct);
	*/

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
			// fprintf( stdout, "%d vec[%d]\t: %lf ( %lf - %lf )\n", i, j, vec[j], neb->images[i + 1][j], neb->images[i - 1][j] );
			absvec += vec[j] * vec[j];
		}
		absvec = sqrt(absvec);
		for(j = 0;j < neb->argdim;j++) { neb->tanmat[i][j] = vec[j] / absvec; }
	}

	free( vec );

	/* chk */
	/*
	fprintf( stdout, "tanmat\n" );
	for(i = 1;i < neb->nimage - 1;i++) { 
		for(j = 0;j < neb->argdim;j++) {
			fprintf( stdout, "%d %d %lf\t", i, j, neb->tanmat[i][j] );
		}
		printf("\n");
	}
	*/
	/* ok */

} // return set of tangent vector


double* Force(NEBInfo* neb, int image) { // printf("Force %d start\n",num);
	int i;
	double *Fvec, *Fs, *Fs_projection, *dV, *dV_projection;

	Fvec = (double*)malloc(sizeof(double) * neb->argdim);
	Fs = Force_spring( neb, image ); // malloc // ok
	Fs_projection = projection(neb->argdim, Fs, neb->tanmat[image] ); // malloc // here
	dV = d_ModelFunc( neb, image ); // malloc // here
	dV_projection = projection(neb->argdim, dV, neb->tanmat[image]); // malloc // ok

	for(i = 0;i < neb->argdim;i++) Fvec[i] = Fs_projection[i] - dV_projection[i]; // hehe NAN 

	/* chk */
	/*
	for(i = 0;i < neb->argdim;i++) {
		fprintf( stdout, "Fs[%d] (%lf), dV[%d] (%lf)\n", i, Fs[i], i, dV[i] );
		fprintf( stdout, "Fvec[%d] (%lf)  = Fs_projection[%d] (%lf) - dV_projection[%d] (%lf)\n", i, Fvec[i], i, Fs_projection[i], i, dV_projection[i]);
	}
	*/

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
	double *dvec = NULL, **dmat = NULL; // grad vector
	FILE *fp;
	char fname[256];

	dmat = (double**)malloc(sizeof(double*) * neb->nimage);
	for(i = 0;i < neb->nimage;i++) { dmat[i] = ( double* )malloc( sizeof( double* ) * neb->argdim ); }

	for(cycle = 0;cycle < neb->cycle;cycle++) { 
		fprintf(stdout, "OptCycle %d start\n", cycle );

		chk = 1; // check converged
		SetTanMat( neb ); // here

		for(i = 1;i < neb->nimage - 1;i++) { // for each images
			dvec = Force(neb,i); // get gradient // malloc
			for(j = 0;j < neb->argdim;j++) dmat[i][j] = dvec[j];
			for(j = 0;j < neb->argdim;j++) { if(dvec[j] > neb->threshold) chk = 0; } // chk converge

			/* chk */
			/*
			fprintf( stdout, "dvec[%d]\t: ", i );
			for(j = 0;j < neb->argdim;j++) { fprintf( stdout, "%lf ", dvec[j] ); }
			fprintf( stdout, "\n" );
			*/

			free(dvec);
		}

		if(chk == 1) break; // converged

		for(i = 1;i < neb->nimage - 1;i++) { 
			for(j = 0;j < neb->argdim;j++) { neb->images[i][j] -= (neb->alpha) * dmat[i][j]; }
		} // update coordinate

		if(cycle %10 == 0) { // FILE OUT start
			sprintf(fname,"%s_NEB_cycle%d.neb.xyz", neb->name, cycle);
			fp = fopen( fname, "w" );
			for(i = 0;i < neb->nimage;i++) {
				fprintf( fp, "%d\n%d\n" ,neb->natom, i);
				for(j = 0;j < neb->natom;j++) { fprintf( fp, "%s\t%17lf\t%17lf\t%17lf\n", neb->elms[j], neb->images[i][3*j], neb->images[i][3*j + 1], neb->images[i][3*j + 2] ); }
			}
			fclose(fp); // FILE OUT end
		} // FILE OUT

	} // Optimization (steepest decent)

	for(i = 0;i < neb->nimage;i++) free( dmat[i] );
	free( dmat );
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

	neb->k = 1.0;
	neb->d = 0.0001;
	neb->threshold = 0.0001;
	neb->alpha = 0.1;

	sprintf( infname, "%s.in", neb->name );
	fp = fopen( infname, "r" );
	if( !fp ) { printf("%s not found\n", infname ); return; }
	while( fgets( line, 256, fp ) ) {
		pt = strstr( line, ":" );
		if( strstr( line, "pes") ) sscanf( pt + 2, "%s", neb->pes );
		else if( strstr( line, "ini") ) sscanf( pt + 2, "%s", neb->xyz );
		else if( strstr( line, "log") ) sscanf( pt + 2, "%s", neb->log );
		else if( strstr( line, "image") ) sscanf( pt + 2, "%d", &neb->nimage );
		else if( strstr( line, "atom") ) sscanf( pt + 2, "%d", &neb->natom );
		else if( strstr( line, "order") ) sscanf( pt + 2, "%d", &neb->order );
		else if( strstr( line, "sconst") ) sscanf( pt + 2, "%lf", &neb->k ); // spring const.
		else if( strstr( line, "ssize") ) sscanf( pt + 2, "%lf", &neb->d ); // step size
		else if( strstr( line, "thold") ) sscanf( pt + 2, "%lf", &neb->threshold ); // threshold of optimization
		else if( strstr( line, "alpha") ) sscanf( pt + 2, "%lf", &neb->alpha ); // 
		else if( strstr( line, "cycle") ) sscanf( pt + 2, "%d", &neb->cycle );
	}
	fclose( fp );

	neb->argdim = neb->natom * 3;
	neb->npair = neb->natom * ( neb->natom - 1 ) / 2;
	neb->nbase = Combination( neb->npair + neb->order, neb->order );
}


void LoadPES( NEBInfo* neb ) {
	int i;
	FILE *fp;
	char line[256];

	fp = fopen( neb->pes, "r" );
	for(i = 0;i < neb->nbase;i++) {
		fgets( line, 256, fp );
		sscanf( line, "%lf", &neb->m->coeff[i] );
		// fprintf( stdout, "coeff[%d]\t: %lf\n", i, neb->m->coeff[i] );
	}
	fclose(fp);
}


void OutNEB( NEBInfo* neb ) {
	FILE *fp_neb;
	int i, j;
	
	fp_neb = fopen( neb->log, "w" );
	for(i = 0;i < neb->nimage;i++) {
		fprintf(fp_neb,"%d\n%d\n", neb->natom, i);
		for(j = 0;j < neb->natom;j++) {
			fprintf(fp_neb, "%s\t%17lf\t%17lf\t%17lf\n", neb->elms[j], neb->images[i][3*j], neb->images[i][3*j + 1], neb->images[i][3*j + 2]);
		}
	}
	fclose( fp_neb );

}


void NEBInfo_malloc( NEBInfo *neb, int chk ) {
	int i, j;
	if( chk > 0 ) { // malloc
		neb->imat = MakeCombination( neb->npair, neb->order );
		neb->images = (double**)malloc( sizeof(double*) * neb->nimage );
		for(i = 0;i < neb->nimage;i++) { neb->images[i] = (double*)malloc(sizeof(double) * neb->argdim); }
		neb->tanmat = (double**)malloc( sizeof(double*) * neb->nimage );
		for(i = 0;i < neb->nimage;i++) { neb->tanmat[i] = (double*)malloc(sizeof(double) * neb->argdim); }
		neb->elms = (char**)malloc(sizeof(char*) * neb->natom);
		for(i = 0;i < neb->natom;i++) { neb->elms[i] = (char*)malloc(sizeof(char) * 10); }

		neb->m = (Model*)malloc(sizeof(Model));
		neb->m->coeff = (double*)malloc(sizeof(double) * neb->nbase);
		neb->m->binfo = (BaseInfo**)malloc(sizeof(BaseInfo*) * neb->nbase);
		neb->m->nbase = neb->nbase;

		for(i = 0;i < neb->nbase;i++) { // for each base function, ...
			neb->m->binfo[i] = (BaseInfo*)malloc(sizeof(BaseInfo));
			neb->m->binfo[i]->expvec = neb->imat[i];
			neb->m->binfo[i]->npair = neb->npair;
			neb->m->binfo[i]->fvec = (Func*)malloc(sizeof(Func) * neb->npair);
			for(j = 0;j < neb->npair;j++) { neb->m->binfo[i]->fvec[j] = Exptype; }
		}
	} else { // free
		for(i = 0;i < neb->nbase;i++) {
			free( neb->m->binfo[i] );
			free( neb->m->binfo[i]->fvec );
		}
		
		free( neb->m->binfo );
		free( neb->m->coeff);
		free( neb->m );

		for(i = 0;i < neb->natom;i++) { free( neb->elms[i] ); }
		free( neb->elms );
		for(i = 0;i < neb->nimage;i++) { free( neb->tanmat[i] ); }
		free( neb->tanmat );
		for(i = 0;i < neb->nimage;i++) { free( neb->images[i] ); }
		free( neb->images );
		free( neb->imat );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
	NEBInfo* neb;

	neb = (NEBInfo*)malloc( sizeof(NEBInfo) ); // resore calc. info.
	sprintf( neb->name, "%s", argv[1] );

	LoadInfile( neb ); // Load info. from input file
	NEBInfo_malloc( neb, 1 ); // Alloc memory
	LoadIniStr( neb ); // Load initial strucutres ( start point and end point )
	LoadPES( neb ); // Load info. about potential ( coefficients )
	MakeIGuess( neb );
	Optimization_ModelFunc( neb );
	OutNEB( neb );
	NEBInfo_malloc( neb, -1 ); // Free memory

	free( neb );

	return 0;
}
