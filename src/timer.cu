//
// Created by leo on 7/8/23.
//

#include "timer.cuh"


Timer::Timer() {
    cudaEventCreate(&event_start);
    cudaEventCreate(&event_stop);
}

Timer::~Timer() {
    cudaEventDestroy(event_start);
    cudaEventDestroy(event_stop);
}

void Timer::start(cudaStream_t stream) {
    cudaEventRecord(event_start, stream);
}

float Timer::elapsed_time(cudaStream_t stream) {
    cudaEventRecord(event_stop, stream);
    cudaEventSynchronize(event_stop);
    float elapsed_time;
    cudaEventElapsedTime(&elapsed_time, event_start, event_stop);
    return elapsed_time;
}

