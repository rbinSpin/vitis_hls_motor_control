#ifndef _PYNQ_HLS_SETTINGS_
#define _PYNQ_HLS_SETTINGS_

#include "ap_axi_sdata.h"
#include "ap_int.h"
#include <inttypes.h>

#define N 4
#define N2 16 // N*N

#define DWIDTH 512
typedef ap_axiu<DWIDTH, 0, 0, 0> axis_t;
typedef float DataType;

const int DataTypeSize = sizeof(DataType) * 8;

typedef union converter {
  DataType d;
  uint32_t i;
} converter_t;

template <typename T> void sda_main(T a[N2], T c[N2]);

extern "C" {
    void Riccati_Solver(hls::stream<axis_t> &, hls::stream<axis_t> &);
}

#endif
