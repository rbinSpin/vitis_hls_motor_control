#include "PYNQ_hls_settings.h"
#include "matrix_ops.h"
#include "hls_stream.h"
#include <math.h>
#include <hls_math.h>   // hls::sqrt / hls::fxp_sqrt / hls::sqrtf

template <typename T> void sda_main(T inputArr[N2], T outputArr[N2])
{
	const int r = 925;
  const int step = 6;

  float As[3][3] = {0};
  float Gs[3][3] = {0};
  float Hs[3][3] = {0};
	
  float Init_Var1[3][3] = {0};
  float Init_Var2[3][3] = {0};
  float Init_Var3[3][3] = {0};
	float Init_Var4[3][3] = {0};
	float tempM1[3][3] = {0};
  float tempM2[3][3] = {0};
	float tempM3[3][3] = {0};
	float As_new[3][3] = {0};
	float GsHs[3][3] = {0};
	float IGsHs[3][3] = {0};
	float As_T_Hs[3][3] = {0};
	float As_T[3][3] = {0};
	float Gs_new_temp[3][3] = {0};
	float Gs_new[3][3] = {0};
	float G_inv_T_Ar[3][3] = {0};
	float inv_Ar_T_Q[3][3] = {0};
	float Hs_new_temp[3][3] = {0};
	float Hs_new[3][3] = {0};
	float Ar[3][3] = {0};
	float Ar_T[3][3] = {0};
	float As_new_temp[3][3] = {0};

	float Var1[3][3] = {0};
	float Var2[3][3] = {0};
	
	float Identity_3d[3][3] = {{1,0,0},
							   {0,1,0},
							   {0,0,1}};
	float Identity_4d[4][4] = {{1,0,0,0},
							   {0,1,0,0},
							   {0,0,1,0},
							   {0,0,0,1}};


	const float Ld = 0.018;             // Inductance of motor [H]
	const float Lq = 0.032;             // Inductance of motor [H]

	float T_e = 0 , omega_hat = 0;
	float omega = 0, iq = 0, id = 0;
	float Vq = 0, Vd = 0;

	float u[2][3] = {0};

	float state_error[3] = {0,0,0};

	float id_hat = 0, iq_hat = 0;


	const float k1 = 5209.23913;		// k1 = 3*p*p/2/J/4*Lambdam ;
	const float k2 = 29.117647;			// k2 = B_value/J ;
	const float k3 = 3836.31713;		// k3 = p/2/J ;
	const float k4 = 78.125;			// k4 = Rs/Lq ;
	const float k5 = 9.4296875;			// k5 = Lambdam/Lq ;
	const float k6 = 31.25;				// k6 = 1/Lq ;
	const float k7 = 138.8888889;		// k7 = Rs/Ld ;
	const float k8 = 55.55555566;		// k8 = 1/Ld ;
	const float k9 = 1.777777778;		// k9 = Lq/Ld ;
	const float k10 = 0.5625;			// k10 = Ld/Lq ;
	const float k11 = -241.68797953;	// k11 = 3*p*p/2/J/4*(Ld-Lq) ;
	const float k12 = -0.0463960;		// k12 = (Ld-Lq)/Lambdam ;
			
//////////////////////	input  //////////////////////
	T_e = inputArr[0];
	omega_hat = inputArr[1];

	omega = inputArr[2];
	iq = inputArr[3];
	id = inputArr[4];

	//Vq, Vd
	Vq = inputArr[5];
	Vd = inputArr[6];

//////////////////////	input  //////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
//                    	MTPA
/////////////////////////////////////////////////////////////////////////////////////////////////////
	const float p = 464.556906;             // p = Lambdam**2/(Ld-Lq)**2 ;
	float q = -342.120181406 * T_e;   // q = ((-4*Lambdam)/(3*poles(Ld-Lq)**2)) * T_e ;
	float delta = hls::powrf(q/2,2) + hls::powrf(p/3,3);

	iq_hat = hls::powrf(-1/q - hls::sqrtf(delta), 1.0f/3.0f) + hls::powrf(-1/q + hls::sqrtf(delta), 1.0f/3.0f);;
	id_hat = k12*iq_hat*iq_hat;

	state_error[0] = omega - omega_hat;
	state_error[1] = iq - iq_hat;
	state_error[2] = id - id_hat;
	
/////////////////////////////////////////////////////////////////////////////////////////////////////
//                    Controller
/////////////////////////////////////////////////////////////////////////////////////////////////////

	#pragma HLS ARRAY_PARTITION variable=state_error type=complete
	
	//    G = B*inv(R)*Transpose(B)
    float G[3][3] = {{0,        0,          0},
					 {0, 976562.5,          0},
					 {0,        0, 3086419.75}};
	
	//	invR*Bc_T
	const float u_temp[2][3] = {{0, -31250.0,           0},
							   	{0,        0, -55555.5556}};
	
	float Q[3][3] = {{1, 0, 0},
					 {0, 5, 0},
					 {0, 0, 5}};      //Qc
	

	// Ar_con = Ar = A-r*I
	Ar[0][0] = -k2 - r;
	Ar[0][1] = k1 + k11*id_hat;
	Ar[0][2] = k11*iq;
	Ar[1][0] = -k5 - k10*id;
	Ar[1][1] = -k4-r;
	Ar[1][2] = -k10*omega_hat;
	Ar[2][0] = k9*iq;
	Ar[2][1] = k9*omega_hat;
	Ar[2][2] = -k7-r;

	Transpose_3d(Ar, Ar_T);

	//////////////////////////// Controller SDA /////////////////////////////////////////////////////////

	// init:
	// inv(Transpose(Ar))*H
	Gaussian_Elimination_3d(Ar_T, Q, inv_Ar_T_Q);

	// G*inv(Transpose(Ar))*H
	Multiply_3d(G, inv_Ar_T_Q, Init_Var1);

	// Ar+G*inv(Transpose(Ar))*H
	Add_3d(Ar, Init_Var1, Init_Var2);

	// inv(Ar+G*inv(Transpose(Ar))*H)
	Gaussian_Elimination_3d(Init_Var2, Identity_3d, tempM1);

	// Transpose(Ar+G*inv(Transpose(Ar))*H)
	Transpose_3d(Init_Var2, Init_Var3);

	// Transpose(inv(Transpose(Ar))*H)
	Transpose_3d(inv_Ar_T_Q, Init_Var4);

	// inv(Transpose(Ar))*H*inv(Ar+G*inv(Transpose(Ar))*H)
	Gaussian_Elimination_T_3d(Init_Var3, Init_Var4, tempM2);

	// G*inv(Transpose(Ar))
	Gaussian_Elimination_T_3d(Ar, G, G_inv_T_Ar);

	// inv(Ar+G*inv(Transpose(Ar))*H)*G*inv(Transpose(Ar))
	Gaussian_Elimination_3d(Init_Var2, G_inv_T_Ar, tempM3);

//	Multiply_3d(tempM1, G_inv_T_Ar, tempM3);
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			As[i][j] = Identity_3d[i][j] + 2*r*tempM1[i][j];
			Hs[i][j] = 2*r*tempM2[i][j];
			Gs[i][j] = 2*r*tempM3[i][j];
		}
	}


	//////////// loop till step //////////////////////////////////////////////////////////////////////
	for(int loop = 0; loop < step; loop++)
	{
		// Gs*Hs
		Multiply_3d(Gs, Hs, GsHs);

		//	I+Gs*Hs
		Add_3d(GsHs, Identity_3d, IGsHs);


		////////// As_new = As*inv(I+Gs*Hs)*As /////////////////////////////////////////////////////////////////////////////
		// inv(I+Gs*Hs)*As
		Gaussian_Elimination_3d(IGsHs, As, As_new_temp);

		// As*inv(I+Gs*Hs)*As
		Multiply_3d(As, As_new_temp, As_new);


		////////// Gs_new = Gs + As*inv(I+Gs*Hs)*Gs*Transpose(As) //////////////////////////////////////////////////////////
		// Transpose(As)
		Transpose_3d(As, As_T);

		// inv(I+Gs*Hs)*Gs
		Gaussian_Elimination_3d(IGsHs, Gs, Var1);

		// inv(I+Gs*Hs)*Gs*Transpose(As)
		Multiply_3d(Var1, As_T, Var2);

		// As*inv(I+Gs*Hs)*Gs*Transpose(As)
		Multiply_3d(As, Var2, Gs_new_temp);

		// Gs + As*Gs*inv(I+Hs*Gs)*Transpose(As)
		Add_3d(Gs, Gs_new_temp, Gs_new);


		////////// Hs_new = Hs + Transpose(As)*Hs*inv(I+Gs*Hs)*As //////////////////////////////////////////////////////////
		// Transpose(As)*Hs
		Multiply_3d(As_T, Hs, As_T_Hs);

		// Transpose(As)*Hs*inv(I+Gs*Hs)
		Multiply_3d(As_T_Hs, As_new_temp, Hs_new_temp);

		// Hs + Transpose(As)*inv(I+Hs*Gs)*Hs*As
		Add_3d(Hs, Hs_new_temp, Hs_new);


		////////// update As,Gs,Hs ///////////////////////////////////////////////////////////////////////

		for(int m = 0; m < 3; m++){
#pragma HLS UNROLL
			for(int n = 0; n < 3; n++){
				As[m][n] = As_new[m][n];
				Gs[m][n] = Gs_new[m][n];
				Hs[m][n] = Hs_new[m][n];
			}
		}
	}
	
	// Hs is the solve of SDA: P


///////////////////////////////////////////////////////////////////////////////////////////////////// 

	for(int i=0; i < 2 ; i++)            //2 times
	{
//#pragma HLS UNROLL
		for(int j=0; j < 3 ; j++)        //3 times
		{
			float temp = 0;
			for(int k=0; k < 3 ; k++)    //3 times
			{
				temp += u_temp[i][k]*Hs[k][j];
			}
			u[i][j] = temp;
		}
	}
//	Multiply_3d(u_temp, Hs, u);

#pragma HLS ARRAY_PARTITION variable=u type=complete
	float Vs[2];

	for(int i=0; i<2; i++)
	{
		float temp = 0;
		for(int j=0;j<3;j++)
		{
			temp += u[i][j]*state_error[j];
		}
		Vs[i] = temp;
	}

	Vq = Vs[0] + (k4*iq_hat + k5*omega_hat + k10*id_hat*omega_hat)*Lq;
	Vd = Vs[1] + (k7*id_hat - k9*omega_hat*iq_hat)*Ld;

	///////////////////////////The above is the Controller///////////////////////////
	///////////////////////////The above is the Controller///////////////////////////
	///////////////////////////The above is the Controller///////////////////////////

	Saturate(Vq, 0.3, -0.3, Vq);
	Saturate(Vd, 0.1, -0.1, Vd);

	outputArr[0] = Vq;
	outputArr[1] = Vd;

	return;
}




extern "C" {
	void Riccati_Solver(hls::stream<axis_t> &in, hls::stream<axis_t> &out) {
//	#pragma HLS INTERFACE s_axilite port = return bundle = control
	#pragma HLS INTERFACE ap_ctrl_none port = return
	#pragma HLS INTERFACE axis port = in
	#pragma HLS INTERFACE axis port = out

	  DataType l_A[N2];
	  DataType l_C[N2];

	#pragma HLS ARRAY_PARTITION variable = l_A factor = 16 dim = 1 cyclic
	#pragma HLS ARRAY_PARTITION variable = l_C factor = 16 dim = 1 cyclic

	  int j_limit = 16; // 512 / DataTypeSize
	  converter_t converter_load[16];
	  converter_t converter_write[16];

	load_A:
		axis_t temp = in.read();
		for(int j = 0; j < j_limit; j++) {
			#pragma HLS UNROLL
			int high = j * DataTypeSize + DataTypeSize - 1;
			int low = j * DataTypeSize;

			converter_load[j].i = temp.data.range(high, low);
			l_A[j] = converter_load[j].d;
		}

		sda_main<DataType>(l_A, l_C);

	writeC:
		axis_t W_temp;
		for(int j = 0; j < j_limit; j++) {
			#pragma HLS UNROLL
			int high = j * DataTypeSize + DataTypeSize - 1;
			int low = j * DataTypeSize;
			converter_write[j].d = l_C[j];
			W_temp.data.range(high, low) = converter_write[j].i;
		}

		W_temp.last = 1;
		W_temp.keep = -1; // enabling all bytes
		out.write(W_temp);
	}
}
