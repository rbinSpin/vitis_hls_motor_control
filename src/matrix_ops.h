#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

template<typename T>
void Swap(T &A, T &B)
{
	T temp = A;
	A = B;
	B = temp;
	return;
}

template<typename T>
void Gaussian_Elimination_3d(T MatrixA[3][3], T MatrixB[3][3], T Minv[3][3])
{
	int j;
	T localVariable, localVariable_copy, localVariable2, localVariable2_copy;
	T MatA[3][3] = {0};
	T MatB[3][3] = {0};

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
        if(MatA[k][k]== (T)0)
        {
        	j = k + 1;
        	while(MatA[j][k] == (T)0)
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


        localVariable = (T)1/MatA[k][k];
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

template<typename T>
void Gaussian_Elimination_4d(T MatrixA[4][4], T MatrixB[4][4], T Minv[4][4])
{
	int j;
	T localVariable, localVariable_copy, localVariable2, localVariable2_copy;
	T MatA[4][4] = {0};
	T MatB[4][4] = {0};

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
        if(MatA[k][k]==(T)0)
        {
        	j = k + 1;
        	while(MatA[j][k] == (T)0)
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


        localVariable = (T)1/MatA[k][k];
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

template<typename T>
void Gaussian_Elimination_T_3d(T MatrixA[3][3], T MatrixB[3][3], T Minv[3][3])
{
	int j;
	T localVariable, localVariable_copy, localVariable2, localVariable2_copy;
	T MatA[3][3] = {0};
	T MatB[3][3] = {0};

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
        if(MatA[k][k]==(T)0)
        {
        	j = k + 1;
        	while(MatA[j][k] == (T)0)
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


        localVariable = (T)1/MatA[k][k];
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

template<typename T>
void Gaussian_Elimination_T_4d(T MatrixA[4][4], T MatrixB[4][4], T Minv[4][4])
{
	int j;
	T localVariable, localVariable_copy, localVariable2, localVariable2_copy;
	T MatA[4][4] = {0};
	T MatB[4][4] = {0};

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
        if(MatA[k][k]==(T)0)
        {
        	j = k + 1;
        	while(MatA[j][k] == (T)0)
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

        localVariable = (T)1/MatA[k][k];
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

template<typename T>
void Multiply_3d(T MatA[3][3], T MatB[3][3], T MatC[3][3])
{
//#pragma HLS ARRAY_PARTITION variable=MatB type=complete
	for (int i = 0; i < 3; i++) {
	    for (int k = 0; k < 3; k++) {
	    	T temp = 0;
	        for (int j = 0; j < 3; j++) {
	        	temp += MatA[i][j]*MatB[j][k];
	        }
	        MatC[i][k] = temp;
	    }
	}
	return;
}

template<typename T>
void Multiply_4d(T MatA[4][4], T MatB[4][4], T MatC[4][4])
{
//#pragma HLS ARRAY_PARTITION variable=MatB type=complete
	for (int i = 0; i < 4; i++) {
	    for (int k = 0; k < 4; k++) {
	    	T temp = 0;
	        for (int j = 0; j < 4; j++) {
	        	temp += MatA[i][j]*MatB[j][k];
	        }
	        MatC[i][k] = temp;
	    }
	}
	return;
}

template<typename T>
void Add_3d(T MatA[3][3], T MatB[3][3], T MatC[3][3])
{
#pragma HLS PIPELINE
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			MatC[i][j] = MatA[i][j] + MatB[i][j];
		}
	}
	return;
}

template<typename T>
void Add_4d(T MatA[4][4], T MatB[4][4], T MatC[4][4])
{
#pragma HLS PIPELINE
	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			MatC[i][j] = MatA[i][j] + MatB[i][j];
		}
	}
	return;
}

template<typename T>
void Transpose_3d(T MatA[3][3], T MatB[3][3])
{
#pragma HLS PIPELINE
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			MatB[j][i] = MatA[i][j];
		}
	}
	return;
}

template<typename T>
void Transpose_4d(T MatA[4][4], T MatB[4][4])
{
#pragma HLS PIPELINE
	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			MatB[j][i] = MatA[i][j];
		}
	}
	return;
}

template<typename T>
void Saturate(T A,  T upper_limit, T lower_limit, T &B)
{
	B = (A>upper_limit)?upper_limit:(A<(lower_limit))?(lower_limit):A;
	return;
}


#endif // MATRIX_OPS_H