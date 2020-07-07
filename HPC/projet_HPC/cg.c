/*
 * Sequential implementation of the Conjugate Gradient Method.
 *
 * Authors : Lilia Ziane Khodja & Charles Bouillaguet
 *
 * v1.02 (2020-04-3)
 *
 * CHANGE LOG:
 *    v1.01 : fix a minor printing bug in load_mm (incorrect CSR matrix size)
 *    v1.02 : use https instead of http in "PRO-TIP"
 *
 * USAGE:
 * 	$ ./cg --matrix bcsstk13.mtx                # loading matrix from file
 *      $ ./cg --matrix bcsstk13.mtx > /dev/null    # ignoring solution
 *	$ ./cg < bcsstk13.mtx > /dev/null           # loading matrix from stdin
 *      $  zcat matrix.mtx.gz | ./cg                # loading gziped matrix from
 *      $ ./cg --matrix bcsstk13.mtx --seed 42      # changing right-hand side
 *      $ ./cg --no-check < bcsstk13.mtx            # no safety check
 *
 * PRO-TIP :
 *      # downloading and uncompressing the matrix on the fly
 *	$ curl --silent https://hpc.fil.cool/matrix/bcsstk13.mtx.gz | zcat | ./cg
 */

#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>

#include "mmio.h"
#include <mpi.h>

#define THRESHOLD 1e-8		// maximum tolerance threshold
#define ROOT 0

struct csr_matrix_t {
	int n;			// dimension
	int nz;			// number of non-zero entries
	int *Ap;		// row pointers
	int *Aj;		// column indices
	double *Ax;		// actual coefficient
};

/*************************** Utility functions ********************************/

/* Seconds (wall-clock time) since an arbitrary point in the past */
double wtime()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1e6;
}

/* Pseudo-random function to initialize b (rumors says it comes from the NSA) */
#define ROR(x, r) ((x >> r) | (x << (64 - r)))
#define ROL(x, r) ((x << r) | (x >> (64 - r)))
#define R(x, y, k) (x = ROR(x, 8), x += y, x ^= k, y = ROL(y, 3), y ^= x)
double PRF(int i, unsigned long long seed)
{
	unsigned long long y = i, x = 0xBaadCafe, b = 0xDeadBeef, a = seed;
	R(x, y, b);
	for (int i = 0; i < 31; i++) {
		R(a, b, i);
		R(x, y, b);
	}
	x += i;
	union { double d; unsigned long long l;	} res;
	res.l = ((x << 2) >> 2) | (((1 << 10) - 1ll) << 52);
	return 2 * (res.d - 1.5);
}

/*************************** Matrix IO ****************************************/

/* Load MatrixMarket sparse symetric matrix from the file descriptor f */
struct csr_matrix_t *load_mm(FILE * f)
{
	MM_typecode matcode;
	int n, m, nnz;

	/* -------- STEP 1 : load the matrix in COOrdinate format */
	double start = wtime();

	/* read the header, check format */
	if (mm_read_banner(f, &matcode) != 0)
		errx(1, "Could not process Matrix Market banner.\n");
	if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode))
		errx(1, "Matrix Market type: [%s] not supported (only sparse matrices are OK)", mm_typecode_to_str(matcode));
	if (!mm_is_symmetric(matcode) || !mm_is_real(matcode))
		errx(1, "Matrix type [%s] not supported (only real symmetric are OK)", mm_typecode_to_str(matcode));
	if (mm_read_mtx_crd_size(f, &n, &m, &nnz) != 0)
		errx(1, "Cannot read matrix size");
	fprintf(stderr, "[IO] Loading [%s] %d x %d with %d nz in triplet format\n", mm_typecode_to_str(matcode), n, n, nnz);
	fprintf(stderr, "     ---> for this, I will allocate %.1f MByte\n", 1e-6 * (40.0 * nnz + 8.0 * n));

	/* Allocate memory for the COOrdinate representation of the matrix (lower-triangle only) */
	int *Ti = malloc(nnz * sizeof(*Ti));
	int *Tj = malloc(nnz * sizeof(*Tj));
	double *Tx = malloc(nnz * sizeof(*Tx));
	if (Ti == NULL || Tj == NULL || Tx == NULL)
		err(1, "Cannot allocate (triplet) sparse matrix");

	/* Parse and load actual entries */
	for (int u = 0; u < nnz; u++) {
		int i, j;
		double x;
		if (3 != fscanf(f, "%d %d %lg\n", &i, &j, &x))
			errx(1, "parse error entry %d\n", u);
		Ti[u] = i - 1;	/* MatrixMarket is 1-based */
		Tj[u] = j - 1;
		/*
		 * Uncomment this to check input (but it slows reading)
		 * if (i < 1 || i > n || j < 1 || j > i)
		 *	errx(2, "invalid entry %d : %d %d\n", u, i, j);
		 */
		Tx[u] = x;
	}

	double stop = wtime();
	fprintf(stderr, "     ---> loaded in %.1fs\n", stop - start);

	/* -------- STEP 2: Convert to CSR (compressed sparse row) representation ----- */
	start = wtime();

	/* allocate CSR matrix */
	struct csr_matrix_t *A = malloc(sizeof(*A));
	if (A == NULL)
		err(1, "malloc failed");
	int *w = malloc((n + 1) * sizeof(*w));
	int *Ap = malloc((n + 1) * sizeof(*Ap));
	int *Aj = malloc(2 * nnz * sizeof(*Ap));
	double *Ax = malloc(2 * nnz * sizeof(*Ax));
	if (w == NULL || Ap == NULL || Aj == NULL || Ax == NULL)
		err(1, "Cannot allocate (CSR) sparse matrix");

	/* the following is essentially a bucket sort */

	/* Count the number of entries in each row */
	for (int i = 0; i < n; i++)
		w[i] = 0;
	for (int u = 0; u < nnz; u++) {
		int i = Ti[u];
		int j = Tj[u];
		w[i]++;
		if (i != j)	/* the file contains only the lower triangular part */
			w[j]++;
	}

	/* Compute row pointers (prefix-sum) */
	int sum = 0;
	for (int i = 0; i < n; i++) {
		Ap[i] = sum;
		sum += w[i];
		w[i] = Ap[i];
	}
	Ap[n] = sum;

	/* Dispatch entries in the right rows */
	for (int u = 0; u < nnz; u++) {
		int i = Ti[u];
		int j = Tj[u];
		double x = Tx[u];
		Aj[w[i]] = j;
		Ax[w[i]] = x;
		w[i]++;
		if (i != j) {	/* off-diagonal entries are duplicated */
			Aj[w[j]] = i;
			Ax[w[j]] = x;
			w[j]++;
		}
	}

	/* release COOrdinate representation */
	free(w);
	free(Ti);
	free(Tj);
	free(Tx);
	stop = wtime();
	fprintf(stderr, "     ---> converted to CSR format in %.1fs\n", stop - start);
	fprintf(stderr, "     ---> CSR matrix size = %.1fMbyte\n", 1e-6 * (24. * nnz + 4. * n));

	A->n = n;
	A->nz = sum;
	A->Ap = Ap;
	A->Aj = Aj;
	A->Ax = Ax;
	return A;
}


/*On prépare la distribution en calculant la taille des paquets pour chaque processus*/
void calcul_repartition(int *ap_nb, int *ap_where, int rank, int nproc, int rows)
{

	if(rank == ROOT)
	{
		int local_rest_rows, local_rows;
		local_rows = ceil((double)rows / nproc);
		local_rest_rows = rows - (nproc-1)*local_rows;

		/*Un processus ne doit pas rien faire ou tout faire !*/
		if (local_rest_rows < 1 || local_rest_rows == rows)
		{
			printf("Error in calcul_repartition !\n");
			exit(7);
		}

		for (int p = 0; p < nproc-1; p++)
		{
			ap_nb[p] = local_rows;
			ap_where[p] = p*local_rows;
		}

		ap_nb[nproc-1] = local_rest_rows;
		ap_where[nproc-1] =(nproc-1)*local_rows;
	}

	MPI_Bcast(ap_nb,nproc,MPI_INT,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ap_where,nproc,MPI_INT,ROOT,MPI_COMM_WORLD);
}


/*Distribution de la matrice*/
void distribution_matrice(const struct csr_matrix_t *A, struct csr_matrix_t *A_local, int* rows,int rank,
	 int nproc, int *ap_nb, int *ap_where)
{

	int nz_local,n_local;
	int *Aj_local,*Ap_local;
	double *Ax_local;

	if(rank==ROOT)
	{
		*rows = A->n;
		n_local = ceil((double)*rows / nproc);
		int ax_nb[nproc], ax_where[nproc];

		for (int p = 0; p < nproc; p++)
		{
			ax_nb[p] = A->Ap[p*n_local + ap_nb[p]] - A->Ap[p * n_local];
			ax_where[p] = A->Ap[p * n_local];
		}

		//ROOT -----> rows
		MPI_Bcast(rows,1,MPI_INT,ROOT,MPI_COMM_WORLD);

		//ROOT -----> nz_local
		MPI_Scatter(ax_nb,1,MPI_INT,&nz_local,1,MPI_INT,ROOT,MPI_COMM_WORLD);

    	//ROOT -----> Ap_local
		Ap_local= (int *)malloc((n_local+1) * sizeof(int));
		MPI_Scatterv(A->Ap,ap_nb,ap_where,MPI_INT,Ap_local,n_local,MPI_INT,ROOT,MPI_COMM_WORLD);
		Ap_local[n_local]=nz_local;

    	//ROOT -----> Ax_local
		Ax_local = (double *)malloc(sizeof(double) * nz_local);
		MPI_Scatterv(A->Ax,ax_nb,ax_where,MPI_DOUBLE,Ax_local,nz_local,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

		//ROOT -----> Aj_local
		Aj_local = (int *)malloc(sizeof(int) * nz_local);
		MPI_Scatterv(A->Aj,ax_nb,ax_where,MPI_INT,Aj_local,nz_local,MPI_INT,ROOT,MPI_COMM_WORLD);

	}else //non ROOT
	{
		//P(i) <----- rows
		MPI_Bcast(rows,1,MPI_INT,ROOT,MPI_COMM_WORLD);

		MPI_Bcast(rows,1,MPI_INT,ROOT,MPI_COMM_WORLD);
		n_local = ceil((double)*rows/ nproc);
		if(rank==nproc-1)
		{
			n_local=*rows - (nproc-1)*n_local;
		}

		//P(i) <----- nz_local
		MPI_Scatter(MPI_IN_PLACE,1,MPI_INT,&nz_local,1,MPI_INT,ROOT,MPI_COMM_WORLD);

		//P(i) <----- Ap_local
		Ap_local= (int *)malloc((n_local+1) * sizeof(int));
		MPI_Scatterv(MPI_IN_PLACE,0,MPI_IN_PLACE,MPI_INT,Ap_local,n_local,MPI_INT,ROOT,MPI_COMM_WORLD);
		Ap_local[n_local] = nz_local;

		//P(i) <----- Ax_local
		Ax_local = (double *)malloc(nz_local * sizeof(double));
		MPI_Scatterv(MPI_IN_PLACE,0,MPI_IN_PLACE,MPI_DOUBLE,Ax_local,nz_local,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

		//P(i) <----- Aj_local
		Aj_local = (int *)malloc(nz_local * sizeof(double));
		MPI_Scatterv(MPI_IN_PLACE,0,MPI_IN_PLACE,MPI_INT,Aj_local,nz_local,MPI_INT,ROOT,MPI_COMM_WORLD);

		//Merci Bouba !! Facilite pour la suite !
		int  prec=Ap_local[0];
		int tmp=0;
		Ap_local[0]=0;
		for(int i=1;i<n_local;i++){
			tmp+=Ap_local[i]-prec;
			prec=Ap_local[i];
			Ap_local[i]=tmp;
		}
	}

	//Initialisation de A_local
	A_local->n = n_local;
	A_local->nz = nz_local;
	A_local->Ap = Ap_local;
	A_local->Aj = Aj_local;
	A_local->Ax = Ax_local;

}


/*************************** Matrix accessors *********************************/
/* Copy the diagonal of A into the vector d. */
void extract_diagonal(const struct csr_matrix_t *A_local, double *d,int rows,int rank,int nproc)
{
	int n_local = A_local->n;
	int *Ap = A_local->Ap;
	int *Aj = A_local->Aj;
	double *Ax = A_local->Ax;

	//Si un proc a un n_local different !
	if( ( (rows % nproc) != 0 ) && rank==nproc-1)
	{
		int start=rows-n_local;
		for (int i = 0; i < n_local; i++) {
			d[i] = 0.0;
			for (int u = Ap[i]; u < Ap[i + 1]; u++){
			if (start+i == Aj[u]){
					d[i] += Ax[u];

				}
			}
		}
	}else
	{
		for (int i = 0; i < n_local; i++) {
				d[i] = 0.0;
				for (int u = Ap[i]; u < Ap[i + 1]; u++)
				if (rank*n_local+i == Aj[u])
						d[i] += Ax[u];
			}
	}

}


/* Matrix-vector product (with A in CSR format) : y = Ax */
void sp_gemv(const struct csr_matrix_t *A, const double *x, double *y)
{
	int n = A->n;
	int *Ap = A->Ap;
	int *Aj = A->Aj;
	double *Ax = A->Ax;

	for (int i = 0; i < n; i++)
	{
		y[i] = 0;
		for (int u = Ap[i]; u < Ap[i + 1]; u++)
		{
			int j = Aj[u];
			double A_ij = Ax[u];
			y[i] += A_ij * x[j];
		}
	}
}


/*************************** Vector operations ********************************/

/* dot product */
double dot(const int n, const double *x, const double *y)
{
	double s = 0.0;
	double res;

	for (int i=0; i < n; i++)
	{
		s += x[i] * y[i];
	}

	MPI_Allreduce(&s, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return res;
}

/* euclidean norm (a.k.a 2-norm) */
double norm(const int n, const double *x)
{
	return sqrt(dot(n, x, x));
}

/*********************** conjugate gradient algorithm *************************/

/* Solve values == b (the solution is written in x). Scratch must be preallocated of size 6n */
void cg_solve(const struct csr_matrix_t *A_local, const double *b_local, double *x, const double epsilon,
			int rows,int rank,int nproc,int* ap_nb,int* ap_where)
{
	int n_local = A_local->n;
	int nz_local = A_local->nz;

  	if(rank==ROOT)
  	{
		fprintf(stderr, "[CG] Starting iterative solver\n");
		fprintf(stderr, "     ---> Working set : %.1fMbyte\n", 1e-6 * (12.0 * nz_local + 52.0 * n_local));
		fprintf(stderr, "     ---> Per iteration: %.2g FLOP in sp_gemv() and %.2g FLOP in the rest\n",
			2. * nz_local, 12. * n_local);
	}

	double *r_local = malloc(n_local * sizeof(double));	      // residue
	double *z_local =malloc(n_local * sizeof(double));	      // preconditioned-residue
	double *p = malloc(rows * sizeof(double));
	double *p_local =malloc(n_local * sizeof(double));	      // search direction
	double *q = malloc(n_local * sizeof(double));	       	  // q == ptr
	double *d_local = malloc(n_local * sizeof(double));	      // diagonal entries of A (Jacobi preconditioning)

	double *x_local = malloc(n_local * sizeof(double));

	/* Isolate diagonal */
	extract_diagonal(A_local, d_local, rows, rank, nproc);

	/*
	 * This function follows closely the pseudo-code given in the (english)
	 * Wikipedia page "Conjugate gradient method". This is the version with
	 * preconditionning.
	 */

	/* We use x_local == 0 --- this avoids the first matrix-vector product. */
	for (int i = 0; i < n_local; i++)
		x_local[i] = 0.0;

	for (int i = 0; i < n_local; i++)	// r <-- b - Ax == b
		r_local[i] = b_local[i];

	for (int i = 0; i < n_local; i++)	// z <-- M^(-1).r
		z_local[i] = r_local[i] / d_local[i];

	for (int i = 0; i < n_local; i++)	// p <-- z
		p_local[i] = z_local[i];


	double rz = dot(n_local, r_local, z_local);
	double start = wtime();
	double last_display = start;
	int iter = 0;

	while (norm(n_local, r_local) > epsilon)
	{
		/* loop invariant : rz = dot(r, z) */
		double old_rz = rz;

		MPI_Allgatherv(p_local, n_local, MPI_DOUBLE, p,ap_nb,ap_where, MPI_DOUBLE, MPI_COMM_WORLD);
		sp_gemv(A_local, p, q);	/* q <-- A.p */

		double alpha = old_rz / dot(n_local, p_local, q);

		for (int i = 0; i < n_local; i++)	// x <-- x + alpha*p
			x_local[i] += alpha * p_local[i];

		for (int i = 0; i < n_local; i++)	// r <-- r - alpha*q
			r_local[i] -= alpha * q[i];

		for (int i = 0; i < n_local; i++)	// z <-- M^(-1).r
			z_local[i] = r_local[i] / d_local[i];

		rz = dot(n_local, r_local, z_local);	// restore invariant
		double beta = rz / old_rz;

		for (int i = 0; i < n_local; i++)		// p <-- z + beta*p
			p_local[i] = z_local[i] + beta * p_local[i];

		iter++;

		if(rank==ROOT)
		{
			double t = wtime();
			if (t - last_display > 0.5) {
				/* verbosity */
				double rate = iter / (t - start);	// iterations per s.
				double GFLOPs = 1e-9 * rate * (2 * nz_local + 12 * n_local);
				fprintf(stderr, "\r     MASTER-->iter : %d (%.1f it/s, %.2f GFLOPs)",iter, rate, GFLOPs);
				fflush(stdout);
				last_display = t;
			}
		}
	}

 	//ROOT <----- P(i) & ROOT
	MPI_Gatherv(x_local,n_local,MPI_DOUBLE,x,ap_nb,ap_where,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

	//On libere
	free(p);free(p_local);free(b_local);free(d_local);free(q);
	free(r_local);free(z_local);free(x_local);


	if(rank==ROOT){
		fprintf(stderr, "\n   ----> Finished in %.1fs and %d  iterations\n", wtime() - start, iter);
	}

}

/******************************* main program *********************************/

/* options descriptor */
struct option longopts[6] = {
	{"seed", required_argument, NULL, 's'},
	{"rhs", required_argument, NULL, 'r'},
	{"matrix", required_argument, NULL, 'm'},
	{"solution", required_argument, NULL, 'o'},
	{"no-check", no_argument, NULL, 'c'},
	{NULL, 0, NULL, 0}
};

int main(int argc, char **argv)
{
	int nproc, rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	struct csr_matrix_t *A=NULL;
	struct csr_matrix_t *A_local=NULL;

	//Utile pour ScatterV plus tard
	int ap_nb[nproc];		//ap_nb[i] <-- la taille du paquet pour P(i)
	int ap_where[nproc];	//ap_where[i] <-- la localisation du paquet pour P(i)
	int n = 0;

	//Tous en ont besoin
	double *b_local;
	double* x={0};
	double* b={0};
	double* scratch=NULL;	/* workspace for cg_solve() */

	/* Parse command-line options */
	int safety_check=0;
	char *solution_filename;


	if(rank==ROOT)
	{
		long long seed = 0;
		char *rhs_filename = NULL;
		char *matrix_filename = NULL;
		solution_filename = NULL;
		safety_check = 1;
		char ch;

		while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
			switch (ch) {
			case 's':
				seed = atoll(optarg);
				break;
			case 'r':
				rhs_filename = optarg;
				break;
			case 'm':
				matrix_filename = optarg;
				break;
			case 'o':
				solution_filename = optarg;
				break;
			case 'c':
				safety_check = 0;
				break;
			default:
				errx(1, "Unknown option");
			}
		}

		/* Load the matrix */
		FILE *f_mat = stdin;
		if (matrix_filename)
		{
			f_mat = fopen(matrix_filename, "r");
			if (f_mat == NULL)
				err(1, "cannot matrix file %s", matrix_filename);
		}
		A = load_mm(f_mat);

		/* Allocate memory */
		n = A->n;
		double *mem = malloc(7 * n * sizeof(double));
		if (mem == NULL)
			err(1, "cannot allocate dense vectors");
		x = mem;				/* solution vector */
		b = mem + n;			/* right-hand side */
		scratch = mem + 2 * n;	/* workspace for cg_solve() */

		/* Prepare right-hand size */
		if (rhs_filename)
		{
			/* load from file */
			FILE *f_b = fopen(rhs_filename, "r");
			if (f_b == NULL)
				err(1, "cannot open %s", rhs_filename);
			fprintf(stderr, "[IO] Loading b from %s\n", rhs_filename);
			for (int i = 0; i < n; i++)
			{
				if (1 != fscanf(f_b, "%lg\n", &b[i]))
					errx(1, "parse error entry %d\n", i);
			}
			fclose(f_b);
		} else
		{
			for (int i = 0; i < n; i++)
				b[i] = PRF(i, seed);
		}

	}//FIN ROOT

	//Calculer la répartition et l'envoyer au processus (Ce n'est pas le broadcast de A)
	calcul_repartition(ap_nb, ap_where, rank, nproc, n);

	//Allouer de l'espace pour matrice locale
	A_local = (struct csr_matrix_t*)malloc(sizeof(struct csr_matrix_t));

  	//Distribution de la matrice (C'est le broadcast de A)
	distribution_matrice(A, A_local, &n, rank, nproc, ap_nb, ap_where);

	//Chaque processus a une partie de b (Il faudrait en discuter car ils need b entier un moment )
	b_local=(double *)malloc(A_local->n * sizeof(double));
	MPI_Scatterv(b,ap_nb,ap_where,MPI_DOUBLE,b_local,A_local->n,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

	//Execution en local
	cg_solve(A_local, b_local, x, THRESHOLD, n, rank, nproc, ap_nb, ap_where);

	if(rank==ROOT)
	{
		/* Dump the solution vector */
	  	FILE *f_x = stdout;
		if (solution_filename != NULL)
		{
			f_x = fopen(solution_filename, "w");
			if (f_x == NULL)
				err(1, "cannot open solution file %s", solution_filename);
			fprintf(stderr, "[IO] writing solution to %s\n", solution_filename);
		}

		/*for (int i = 0; i < n; i++)
			fprintf(f_x, "%a\n", x[i]);*/
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}