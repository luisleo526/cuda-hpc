//
// Created by leo on 7/8/23.
//

#include "constants.cuh"

#ifndef CUHPC_TIMER_CUH
#define CUHPC_TIMER_CUH

class Timer {

public:
    cudaEvent_t event_start, event_stop;
    Timer();
    ~Timer();

    void start(cudaStream_t stream = 0);
    float elapsed_time(cudaStream_t stream = 0);


};


#endif //CUHPC_TIMER_CUH
