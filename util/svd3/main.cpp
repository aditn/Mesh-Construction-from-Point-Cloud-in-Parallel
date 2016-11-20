#include <stdio.h>
#include <string.h>
#include "svd3.h"

int main(int argc, char **argv) {
	int e;
	double s[3], A[3][3], U[3][3], S[3][3], V[3][3], USVT[3][3], US[3][3], 
		UTU[3][3], VTV[3][3];
	memset(S,0,sizeof(S));

	do {
		puts("\nEnter 3x3 matrix by row (q to quit):");
		e = readmat3(A);

		putchar('\n');

		if (e == 9) {
			printmat3(A);

			puts("\nComputing SVD ...\n");

			svd3((double *)U, s, (double *)V, (double *)A);

			S[0][0] = s[0];
			S[1][1] = s[1];
			S[2][2] = s[2];

			printmat3(U);
			putchar('\n');
			printmat3(S);
			putchar('\n');
			printmat3(V);
			putchar('\n');

			ata3((double *)UTU, (double *)U);
			printmat3(UTU);
			putchar('\n');

			ata3((double *)VTV, (double *)V);
			printmat3(VTV);
			putchar('\n');

			matmul3((double *)US, (double *)U, (double *)S);
			trans3((double *)V);
			matmul3((double *)USVT, (double *)US, (double *)V);

			printmat3(USVT);
		}
		else if (e > 0) {
			puts("Input error. Try again.");
		}

		fpurge(stdin);

	} while (e > 0);

	return 0;
}

