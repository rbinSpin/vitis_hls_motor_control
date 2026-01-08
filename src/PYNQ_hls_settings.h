#ifndef _PYNQ_HLS_SETTINGS_
#define _PYNQ_HLS_SETTINGS_

#include "ap_int.h" // For ap_fixed
#include <inttypes.h> // For any potential integer types needed

#define N 4
#define N2 16 // N*N

typedef ap_fixed<34, 14, AP_TRN, AP_SAT> DataType;

const int DataTypeSize = 34;

void sda_main(DataType a[N2], DataType c[N2]);

#endif