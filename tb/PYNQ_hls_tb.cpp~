#include <iostream>
#include "hls_stream.h"
#include "ap_axi_sdata.h"
#include "ap_int.h"
#include "PYNQ_hls_settings.h"

#define DWIDTH 512
#define DataType float
#define DataTypeSize 32
#define NUM_ELEMS 16  // 512 / 32 = 16

typedef ap_axiu<DWIDTH, 0, 0, 0> axis_t;


int main() {
    hls::stream<axis_t> in_stream("input_stream");
    hls::stream<axis_t> out_stream("output_stream");

    // Step 1: 準備輸入資料
    DataType input_array[NUM_ELEMS];
    for (int i = 0; i < NUM_ELEMS; i++) {
        input_array[i] = (DataType)(i + 1);  // ex: 1.0, 2.0, ..., 16.0
    }

    // Step 2: 打包成 axis_t 並寫入 stream
    axis_t input_packet;
    ap_uint<DWIDTH> data_word = 0;
    converter_t conv;

    for (int i = 0; i < NUM_ELEMS; i++) {
        conv.d = input_array[i];
        int high = i * DataTypeSize + DataTypeSize - 1;
        int low = i * DataTypeSize;
        data_word.range(high, low) = conv.i;
    }

    input_packet.data = data_word;
    input_packet.last = 1;  // only one beat
    input_packet.keep = -1; // all bytes valid
    in_stream.write(input_packet);

    // Step 3: 呼叫待測函數
    Riccati_Solver(in_stream, out_stream);

    // Step 4: 從 output stream 取回結果
    if (!out_stream.empty()) {
        axis_t output_packet = out_stream.read();
        ap_uint<DWIDTH> out_data = output_packet.data;

        DataType output_array[NUM_ELEMS];
        for (int i = 0; i < NUM_ELEMS; i++) {
            int high = i * DataTypeSize + DataTypeSize - 1;
            int low = i * DataTypeSize;
            conv.i = out_data.range(high, low);
            output_array[i] = conv.d;
        }

        // Step 5: 印出輸出結果
        std::cout << "輸出結果: ";
        for (int i = 0; i < NUM_ELEMS; i++) {
            std::cout << output_array[i] << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "Output stream is empty!" << std::endl;
    }

    return 0;
}
