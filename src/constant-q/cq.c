#include <stdlib.h>
#include <math.h>
#include <complex.h>
 
#include "cxsparse/Include/cs.h"
#include <fftw3.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

 
float* cq_hanning_window(int length) {
    float* result = fftwf_alloc_real(length);
 
    for (int j = 0; j < length; ++j) {
        float elem = 0.5f * (1-(cos(j*2*M_PI / (float)(length - 1))));
        result[j] = elem / (float)length;
    }
 
    return result;
}
 
fftwf_complex* cq_temp_kernel(float Q, int length) {
    fftwf_complex* result = fftwf_alloc_complex(length);
 
    for (int j = 0; j < length; ++j) {
        fftwf_complex elem = cexpf(2*M_PI*_Complex_I*Q*j / (float)length);
        result[j] = elem;
    }
 
    return result;
}
 
void cq_multiply_vector_elem(const float* in, const fftwf_complex* in2,
        fftwf_complex* out, int length) {
    for (int j = 0; j < length; ++j) {
        fftwf_complex elem = in[j] * in2[j];
        out[j] = elem;
    }
}
 
float cq_q(int bins) {
    return 1 / (exp2(1 / (float)bins) - 1);
}
 
int cq_k(int min_freq, int max_freq, int bins) {
    return rint(ceil(bins * log2((float)max_freq / (float)min_freq)));
}
 
int cq_fft_len(float q, int min_freq, int sample_rate) {
    return rint(exp2(ceil(log2(q * sample_rate / (float) min_freq))));
}
 
void cq_swap_matrix_column(float _Complex* matrix, int height, int width,
        const float _Complex* vector, int index) {
    for (int j = 0; j < height; ++j) {
        matrix[j * width + index] = vector[j];
    }
}
 
cs_ci* cq_make_kernel(int min_freq, int max_freq, int sample_rate,
        int bins, int* height, int* width) {
    float Q = cq_q(bins);
    int K = cq_k(min_freq, max_freq, bins);
    int fft = cq_fft_len(Q, min_freq, sample_rate);
 
    fftwf_complex* temp = fftwf_alloc_complex(fft);
    for (int j = 0; j < fft; ++j)
        temp[j] = 0.f + 0.f*_Complex_I;
    fftwf_complex* spec = fftwf_alloc_complex(fft);
    cs_ci* result = cs_ci_spalloc(fft, K, (int)(fft*K*0.01f), 1, 1);
    fftwf_plan plan = fftwf_plan_dft_1d(fft, temp, spec, FFTW_FORWARD,
            FFTW_ESTIMATE);
 
    for (int j = K-1; j >= 0; --j) {
        int len = ceil(Q*sample_rate/(min_freq*exp2(j/(float)bins)));
 
        float* h = cq_hanning_window(len);
        fftwf_complex* vec = cq_temp_kernel(Q, len);
 
        // temp = h .* vec
        cq_multiply_vector_elem(h, vec, temp, len);
 
        // spec = fft(temp)
        fftwf_execute(plan);
 
        for (int k = 0; k < fft; ++k) {
            fftwf_complex temp = spec[k];
            if (cabsf(temp) > 0.0054f) {
                cs_ci_entry(result, k, j, conjf(temp) / (float)fft);
            }
        }
 
        fftwf_free(h);
        fftwf_free(vec);
    }
 
    fftwf_destroy_plan(plan);
    fftwf_free(temp);
    fftwf_free(spec);
 
    cs_ci* ret = cs_ci_compress(result);
    cs_ci_spfree(result);
    *height = fft;
    *width = K;
    return ret;
}
 
fftwf_complex* cq_const_q_transform(float* data, const cs_ci* kernel,
        int height, int width, const fftwf_plan plan, fftwf_complex* buffer) {
    fftwf_execute_dft_r2c(plan, data, buffer);
 
    // fill redundant data
    for (int j = 0; j < height / 2 - 1; ++j) {
        buffer[height - j - 1] = conjf(buffer[j + 1]);
    }
 
    fftwf_complex* result = fftwf_alloc_complex(width);
 
    for (int j = 0; j < width; ++j) {
        result[j] = 0.0f;
    }
 
    cs_ci_gaxpy(kernel, buffer, result);
 
    return result;
}
 
fftwf_complex* cq_const_q_wrap(float* data, const cs_ci* kernel, int height,
        int width, const int* indices, int indices_size) {
    fftwf_complex* temp_fft = fftwf_alloc_complex(height);
    fftwf_plan plan = fftwf_plan_dft_r2c_1d(height, data, temp_fft, FFTW_ESTIMATE |
            FFTW_PRESERVE_INPUT);
 
    fftwf_complex* result = fftwf_alloc_complex(height * indices_size);
 
    for (int j = 0; j < indices_size; ++j) {
        fftwf_complex* cq = cq_const_q_transform(data + indices[j], kernel,
            height, width, plan, temp_fft);
 
        cq_swap_matrix_column(result, width, indices_size, cq, j);
        fftwf_free(cq);
    }
 
    fftwf_free(temp_fft);
    fftwf_destroy_plan(plan);
    return result;
}
 
fftwf_complex* cq_short_time_constq_transform(float* data, int data_length,
        int min_freq, int max_freq, int sample_rate, int bins, int step,
        int* height,  int* width) {
    int kernel_height, kernel_width;
    cs_ci* ker = cq_make_kernel(min_freq, max_freq, sample_rate, bins,
        &kernel_height, &kernel_width);
    cs_ci* kern = cs_ci_transpose(ker, 1);
 
    int max_index = rint(ceil(data_length / (float)kernel_height));
 
    int indices_size = (max_index - 1) * kernel_height / (double) step + 1;
    int* indices = (int*) malloc(indices_size * sizeof(int));
    for (int j = 0, k = 0; j <= (max_index - 1) * kernel_height; ++k, j += step)
        indices[k] = j;
 
    fftwf_complex* result = cq_const_q_wrap(data, kern, kernel_height,
            kernel_width, indices, indices_size);
 
    cs_ci_spfree(ker);
    cs_ci_spfree(kern);
    free(indices);
 
    *height = kernel_width;
    *width = indices_size;
    return result;
}