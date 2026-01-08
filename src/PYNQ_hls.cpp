#include "PYNQ_hls_settings.h"
#include "matrix_ops.h"


#include <hls_math.h>   // hls::sqrt / hls::fxp_sqrt / hls::sqrtf
#include "ap_fixed.h"
#include "hls_math.h"

void sda_main(DataType inputArr[N2], DataType outputArr[N2]){
	const int r = 925;
  const int step = 6;

  DataType As[3][3] = {0};
  DataType Gs[3][3] = {0};
  DataType Hs[3][3] = {0};
	
  DataType Init_Var1[3][3] = {0};
  DataType Init_Var2[3][3] = {0};
  DataType Init_Var3[3][3] = {0};
	DataType Init_Var4[3][3] = {0};
	DataType tempM1[3][3] = {0};
  DataType tempM2[3][3] = {0};
	DataType tempM3[3][3] = {0};
	DataType As_new[3][3] = {0};
	DataType GsHs[3][3] = {0};
	DataType IGsHs[3][3] = {0};
	DataType As_T_Hs[3][3] = {0};
	DataType As_T[3][3] = {0};
	DataType Gs_new_temp[3][3] = {0};
	DataType Gs_new[3][3] = {0};
	DataType G_inv_T_Ar[3][3] = {0};
	DataType inv_Ar_T_Q[3][3] = {0};
	DataType Hs_new_temp[3][3] = {0};
	DataType Hs_new[3][3] = {0};
	DataType Ar[3][3] = {0};
	DataType Ar_T[3][3] = {0};
	DataType As_new_temp[3][3] = {0};

	DataType Var1[3][3] = {0};
	DataType Var2[3][3] = {0};
	
	DataType Identity_3d[3][3] = {{ (DataType)1, (DataType)0, (DataType)0},
							   { (DataType)0, (DataType)1, (DataType)0},
							   { (DataType)0, (DataType)0, (DataType)1}};
	DataType Identity_4d[4][4] = {{ (DataType)1, (DataType)0, (DataType)0, (DataType)0},
							   { (DataType)0, (DataType)1, (DataType)0, (DataType)0},
							   { (DataType)0, (DataType)0, (DataType)1, (DataType)0},
							   { (DataType)0, (DataType)0, (DataType)0, (DataType)1}};


	const DataType Ld = 0.018;             // Inductance of motor [H]
	const DataType Lq = 0.032;             // Inductance of motor [H]

	DataType T_e = 0 , omega_hat = 0;
	DataType omega = 0, iq = 0, id = 0;
	DataType Vq = 0, Vd = 0;

	DataType u[2][3] = {0};

	DataType state_error[3] = {0,0,0};

	DataType id_hat = 0, iq_hat = 0;


	const DataType k1  = (DataType)0.8295;       // k1  = (3*p*lambda_m)/(2*J*wb) * (Ib/Ib) normalized
	const DataType k2  = (DataType)2.9118;       // k2  = B/J (Mechanical damping ratio)
	const DataType k3  = (DataType)0.6109;       // k3  = (p*Ib)/(2*J*wb) torque-to-speed coefficient
	const DataType k4  = (DataType)78.1250;      // k4  = R/Lq (Electrical time constant Q)
	const DataType k5  = (DataType)592.1844;     // k5  = (lambda_m*wb)/Vb / (Lq*Ib/Vb) back-EMF term
	const DataType k6  = (DataType)558.8125;     // k6  = Vb/(Lq*Ib) reciprocal of normalized Lq
	const DataType k7  = (DataType)138.8889;     // k7  = R/Ld (Electrical time constant D)
	const DataType k8  = (DataType)993.4444;     // k8  = Vb/(Ld*Ib) reciprocal of normalized Ld
	const DataType k9  = (DataType)1116.4;       // k9  = Lq/Ld ratio scaled by base values
	const DataType k10 = (DataType)353.2500;     // k10 = Ld/Lq ratio scaled by base values
	const DataType k11 = (DataType)-3.8485;      // k11 = Reluctance torque coefficient
	const DataType k12 = (DataType)-0.0463960;		// k12 = (Ld-Lq)/Lambdam ;
			
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
	const DataType p = (DataType)464.556906;             // p = Lambdam**2/(Ld-Lq)**2 ;
	DataType q = (DataType)-342.120181406 * T_e;   // q = ((-4*Lambdam)/(3*poles(Ld-Lq)**2)) * T_e ;
	DataType temp_q = q/(DataType)2;
	DataType temp_p = p/(DataType)3;
	DataType delta = (temp_q*temp_q) + (temp_p*temp_p*temp_p);


	DataType sqrt_delta = hls::sqrtf((DataType)delta);
	iq_hat = hls::cbrt( (DataType)-q/2 - sqrt_delta) + hls::cbrt((DataType)-q/2 + sqrt_delta);
	id_hat = k12*iq_hat*iq_hat;

	state_error[0] = omega - omega_hat;
	state_error[1] = iq - iq_hat;
	state_error[2] = id - id_hat;
	
/////////////////////////////////////////////////////////////////////////////////////////////////////
//                    Controller
/////////////////////////////////////////////////////////////////////////////////////////////////////

	#pragma HLS ARRAY_PARTITION variable=state_error type=complete
	
	//    G = B*inv(R)*B'
	//    B =   [       0         0
	//	 		 558.8125         0
	//					0  993.4444]
	//    R = diag([30, 30])
    DataType G[3][3] = {{(DataType)0,        (DataType)0,          (DataType)0},
					 {(DataType)0, 			(DataType)1041,          (DataType)0},
					 {(DataType)0,       	 (DataType)0, 			(DataType)3290}};
	
	//	invR*B'
	const DataType u_temp[2][3] = {{(DataType)0, (DataType)-31250.0,           (DataType)0},
							   	{(DataType)0,        (DataType)0, (DataType)-55555.5556}};
	
	DataType Q[3][3] = {{(DataType)5000, (DataType)0, (DataType)0},
					 {(DataType)0, (DataType)10, (DataType)0},
					 {(DataType)0, (DataType)0, (DataType)10}};      //Qc
	

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
				DataType temp = 0;
			for(int k=0; k < 3 ; k++)    //3 times
			{
				temp += u_temp[i][k]*Hs[k][j];
			}
			u[i][j] = temp;
		}
	}
//	Multiply_3d(u_temp, Hs, u);

#pragma HLS ARRAY_PARTITION variable=u type=complete
	DataType Vs[2];

	for(int i=0; i<2; i++)
	{
			DataType temp = 0;
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

	Saturate(Vq, (DataType)1.1, (DataType)-1.1, Vq);
	Saturate(Vd, (DataType)1.1, (DataType)-1.1, Vd);

	outputArr[0] = Vq;
	outputArr[1] = Vd;

	return;

}

