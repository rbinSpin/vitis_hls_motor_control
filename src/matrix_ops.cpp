#include "matrix_ops.h"


void Gaussian_Elimination_3d(float MatrixA[3][3], float MatrixB[3][3], float Minv[3][3])
{
	int j;
	float localVariable, localVariable_copy, localVariable2, localVariable2_copy;
	float MatA[3][3] = {0};
	float MatB[3][3] = {0};

	for(int m=0; m<3; m++)
	{
		for(int n=0; n<3; n++)
		{
			MatA[m][n] = MatrixA[m][n];
			MatB[m][n] = MatrixB[m][n];
		}
	}

	for(int k=0; k<3; k++)
	{
		// if diagonal is zero, then return if all row elements are 0
        if(MatA[k][k]==0)
        {
        	j = k + 1;
        	while(MatA[j][k] == 0)
        	{
                if(j >= 2)
                {
    				return;
                }
                else
                {
            		j = j + 1;
                }
			}

            for(int x=k; x < 3; x++)
            {
				#pragma HLS UNROLL
            	Swap(MatA[k][x], MatA[j][x]);
			}

            for(int x=0; x < 3; x++)
            {
				#pragma HLS UNROLL
            	Swap(MatB[k][x], MatB[j][x]);
			}
        }


        localVariable = 1/MatA[k][k];
		localVariable_copy = localVariable;
		MatA[k][k] = 1;

		for(int m=k+1; m<3; m++)
		{
			#pragma HLS UNROLL
			MatA[k][m]*=localVariable;
		}

		for(int m=0; m<3; m++)
		{
			#pragma HLS UNROLL
			MatB[k][m]*=localVariable_copy;
		}


	 	for(int m=0; m<3; m++)
	 	{
	 		if(m != k)
			{
				localVariable2 = MatA[m][k];
				localVariable2_copy = localVariable2;

				MatA[m][k] = 0;

				for(int n=k+1; n<3; n++)
				{
					#pragma HLS UNROLL
					MatA[m][n] -= MatA[k][n]*localVariable2;
				}

				for(int n=0; n<3; n++)
				{
					#pragma HLS UNROLL
					MatB[m][n] -= MatB[k][n]*localVariable2_copy;
				}
			}
	 	}
	}

    for(int m=0; m<3; m++)
    {
		for(int n=0; n<3; n++)
		{
			Minv[m][n] = MatB[m][n];
		}
    }

	return;
}

void Gaussian_Elimination_4d(float MatrixA[4][4], float MatrixB[4][4], float Minv[4][4])
{
	int j;
	float localVariable, localVariable_copy, localVariable2, localVariable2_copy;
	float MatA[4][4] = {0};
	float MatB[4][4] = {0};

	for(int m=0; m<4; m++)
	{
		for(int n=0; n<4; n++)
		{
			MatA[m][n] = MatrixA[m][n];
			MatB[m][n] = MatrixB[m][n];
		}
	}

	for(int k=0;k<4;k++)
	{
		// if diagonal is zero, then return if all row elements are 0
        if(MatA[k][k]==0)
        {
        	j = k + 1;
        	while(MatA[j][k] == 0)
        	{
                if(j >= 3)
                {
    				return;
                }
                else
                {
            		j = j + 1;
                }
			}

            for(int x=k; x < 4; x++)
            {
				#pragma HLS UNROLL
            	Swap(MatA[k][x], MatA[j][x]);
			}

            for(int x=0; x < 4; x++)
            {
				#pragma HLS UNROLL
            	Swap(MatB[k][x], MatB[j][x]);
			}
        }


        localVariable = 1/MatA[k][k];
		localVariable_copy = localVariable;
		MatA[k][k] = 1;

		for(int m=k+1; m<4; m++)
		{
			#pragma HLS UNROLL
			MatA[k][m]*=localVariable;
		}

		for(int m=0; m<4; m++)
		{
			#pragma HLS UNROLL
			MatB[k][m]*=localVariable_copy;
		}


	 	for(int m=0; m<4; m++)
	 	{
	 		if(m != k)
	 		{
				localVariable2 = MatA[m][k];
				localVariable2_copy = localVariable2;

				MatA[m][k] = 0;

				for(int n=k+1; n<4; n++)
				{
					#pragma HLS UNROLL
					MatA[m][n] -= MatA[k][n]*localVariable2;
				}

				for(int n=0; n<4; n++)
				{
					#pragma HLS UNROLL
					MatB[m][n] -= MatB[k][n]*localVariable2_copy;
				}
	 		}
	 	}
	}

    for(int m=0; m<4; m++)
    {
		for(int n=0; n<4; n++)
		{
			Minv[m][n] = MatB[m][n];
		}
    }

	return;
}


void Gaussian_Elimination_T_3d(float MatrixA[3][3], float MatrixB[3][3], float Minv[3][3])
{
	int j;
	float localVariable, localVariable_copy, localVariable2, localVariable2_copy;
	float MatA[3][3] = {0};
	float MatB[3][3] = {0};

	for(int m=0; m<3; m++)
	{
		for(int n=0; n<3; n++)
		{
			MatA[m][n] = MatrixA[m][n];
			MatB[m][n] = MatrixB[m][n];
		}
	}

	for(int k=0;k<3;k++)
	{
		// if diagonal is zero, then return if all row elements are 0
        if(MatA[k][k]==0)
        {
        	j = k + 1;
        	while(MatA[j][k] == 0)
        	{
                if(j >= 2)
                {
    				return;
                }
                else
                {
            		j = j + 1;
                }
			}

            for(int x=k; x < 3; x++)
            {
				#pragma HLS UNROLL
            	Swap(MatA[k][x], MatA[j][x]);
			}

            for(int x=0; x < 3; x++)
            {
				#pragma HLS UNROLL
            	Swap(MatB[k][x], MatB[j][x]);
			}
        }


        localVariable = 1/MatA[k][k];
		localVariable_copy = localVariable;
		MatA[k][k] = 1;

		for(int m=k+1; m<3; m++)
		{
			#pragma HLS UNROLL
			MatA[k][m]*=localVariable;
		}

		for(int m=0; m<3; m++)
		{
			#pragma HLS UNROLL
			MatB[k][m]*=localVariable_copy;
		}


	 	for(int m=0; m<3; m++)
	 	{
	 		if(m != k)
			{
				localVariable2 = MatA[m][k];
				localVariable2_copy = localVariable2;

				MatA[m][k] = 0;

				for(int n=k+1; n<3; n++)
				{
					#pragma HLS UNROLL
					MatA[m][n] -= MatA[k][n]*localVariable2;
				}

				for(int n=0; n<3; n++)
				{
					#pragma HLS UNROLL
					MatB[m][n] -= MatB[k][n]*localVariable2_copy;
				}
			}
	 	}
	}

    for(int m=0; m<3; m++)
    {
		for(int n=0; n<3; n++)
		{
			Minv[m][n] = MatB[n][m];
		}
    }

	return;
}

void Gaussian_Elimination_T_4d(float MatrixA[4][4], float MatrixB[4][4], float Minv[4][4])
{
	int j;
	float localVariable, localVariable_copy, localVariable2, localVariable2_copy;
	float MatA[4][4] = {0};
	float MatB[4][4] = {0};

	for(int m=0; m<4; m++)
	{
		for(int n=0; n<4; n++)
		{
			MatA[m][n] = MatrixA[m][n];
			MatB[m][n] = MatrixB[m][n];
		}
	}

	for(int k=0;k<4;k++)
	{
		// if diagonal is zero, then return if all row elements are 0
        if(MatA[k][k]==0)
        {
        	j = k + 1;
        	while(MatA[j][k] == 0)
        	{
                if(j >= 3)
                {
    				return;
                }
                else
                {
            		j = j + 1;
                }
			}

            for(int x=k; x < 4; x++)
            {
				#pragma HLS UNROLL
            	Swap(MatA[k][x], MatA[j][x]);
			}

            for(int x=0; x < 4; x++)
            {
				#pragma HLS UNROLL
            	Swap(MatB[k][x], MatB[j][x]);
			}
        }

        localVariable = 1/MatA[k][k];
		localVariable_copy = localVariable;
		MatA[k][k] = 1;

		for(int m=k+1; m<4; m++)
		{
			#pragma HLS UNROLL
			MatA[k][m]*=localVariable;
		}

		for(int m=0; m<4; m++)
		{
			#pragma HLS UNROLL
			MatB[k][m]*=localVariable_copy;
		}


	 	for(int m=0; m<4; m++)
	 	{
	 		if(m != k)
	 		{
				localVariable2 = MatA[m][k];
				localVariable2_copy = localVariable2;

				MatA[m][k] = 0;

				for(int n=k+1; n<4; n++)
				{
					#pragma HLS UNROLL
					MatA[m][n] -= MatA[k][n]*localVariable2;
				}

				for(int n=0; n<4; n++)
				{
					#pragma HLS UNROLL
					MatB[m][n] -= MatB[k][n]*localVariable2_copy;
				}
	 		}
	 	}
	}

	for(int m=0; m<4; m++)
	{
		for(int n=0; n<4; n++)
		{
			Minv[m][n] = MatB[n][m];
		}
	}

	return;
}


void Multiply_3d(float MatA[3][3], float MatB[3][3], float MatC[3][3])
{
//#pragma HLS ARRAY_PARTITION variable=MatB type=complete
	for (int i = 0; i < 3; i++) {
	    for (int k = 0; k < 3; k++) {
	    	float temp = 0;
	        for (int j = 0; j < 3; j++) {
	        	temp += MatA[i][j]*MatB[j][k];
	        }
	        MatC[i][k] = temp;
	    }
	}
	return;
}

void Multiply_4d(float MatA[4][4], float MatB[4][4], float MatC[4][4])
{
//#pragma HLS ARRAY_PARTITION variable=MatB type=complete
	for (int i = 0; i < 4; i++) {
	    for (int k = 0; k < 4; k++) {
	    	float temp = 0;
	        for (int j = 0; j < 4; j++) {
	        	temp += MatA[i][j]*MatB[j][k];
	        }
	        MatC[i][k] = temp;
	    }
	}
	return;
}

void Add_3d(float MatA[3][3], float MatB[3][3], float MatC[3][3])
{
#pragma HLS PIPELINE
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			MatC[i][j] = MatA[i][j] + MatB[i][j];
		}
	}
	return;
}

void Add_4d(float MatA[4][4], float MatB[4][4], float MatC[4][4])
{
#pragma HLS PIPELINE
	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			MatC[i][j] = MatA[i][j] + MatB[i][j];
		}
	}
	return;
}

void Transpose_3d(float MatA[3][3], float MatB[3][3])
{
#pragma HLS PIPELINE
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			MatB[j][i] = MatA[i][j];
		}
	}
	return;
}

void Transpose_4d(float MatA[4][4], float MatB[4][4])
{
#pragma HLS PIPELINE
	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			MatB[j][i] = MatA[i][j];
		}
	}
	return;
}

void Swap(float &A,float &B)
{
	float temp = A;
	A = B;
	B = temp;
	return;
}

void Saturate(float A,  float upper_limit, float lower_limit, float &B)
{
	B = (A>upper_limit)?upper_limit:(A<(lower_limit))?(lower_limit):A;
	return;
}
