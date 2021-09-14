#ifndef CQ_FFT_H
#define CQ_FFT_H

#include <complex.h>
#include <fftw3.h>

#include "cxsparse/Include/cs.h"

/**
 *  \brief Prints given vector on stdout
 *
 *  Each value is printed up to 7 decimal places
 *
 *  \param vec Vector to be printed
 *  \param length Length of a given vector
 */
void cq_dump_vector(const double* vec, int length);

/**
 *  \brief Prints given complex vector on stdout
 *
 *  Prints "Re(val) + Im(val)i", each up to 7 decimal places
 *
 *  \param vec Complex vector to be printed
 *  \param length Length of a given vector
 */
void cq_dump_vector_complex(const double _Complex* vec, int length);

/**
 *  \brief Prints given matrix to stdout
 *
 *  Prints magnitude of each value in a matrix
 *
 *  \param matrix Matrix to be printed
 *  \param height
 *  \param width
 */
void cq_dump_matrix_abs(double _Complex* matrix, int height, int width);

/**
 *  \brief Returns hanning window of a given length
 *
 *  Returned vector is normalized by dividing each element by its length
 *
 *  \param length length of a returned window
 */
double* cq_hanning_window(int length);

/**
 *  \brief Returns temporary vector used in computing a transform kernel
 *
 *  \param Q Q value in transform
 *  \param length
 */
fftw_complex* cq_temp_kernel(double Q, int length);

/**
 *  \brief Helper function to compute x = a .* b, where x, a, b are vectors
 *
 *  \param in
 *  \param in2
 *  \param out Vector to put result in, have to be allocated
 *  \param length Length of each vector
 */
void cq_multiply_vector_elem(const double* in, const fftw_complex* in2, fftw_complex* out,
                             int length);

/**
 *  \brief Computes Q coefficient
 *
 *  \param bins Number of bins per octave
 */
double cq_q(int bins);

/**
 *  \brief Computes K coefficient
 *
 *  \param min_freq Lower bound of frequency range
 *  \param max_freq Higher bound of frequency range
 *  \param bins Number of bins per octave
 */
int cq_k(int min_freq, int max_freq, int bins);

/**
 *  \brief Computes length of FFT used with given transform parameters
 *
 *  \param q Q coefficient returned by cq_q(int)
 *  \param min_freq Lower bound of frequency range
 *  \param sample_rate Sample rate of source signal
 */
int cq_fft_len(double q, int min_freq, int sample_rate);

/**
 *  \brief Zeroes values which amplitude is smaller than threshold
 *
 *  \param vec Input vector
 *  \param length Length of a vector
 *  \param thresh Threshold
 */
void cq_zero_vector_below_thresh(double _Complex* vec, int length, double thresh);

/**
 *  \brief Sets values in a given column of a given matrix to those in a vector
 *
 *  \param matrix
 *  \param height Height of a matrix
 *  \param width Width of a matrix
 *  \param vector Vector of values
 *  \param index Index of a column in matrix
 */
void cq_swap_matrix_column(double _Complex* matrix, int height, int width,
                           const double _Complex* vector, int index);

/**
 *  \brief Returns a transform kernel used throughout the algorithm
 *
 *  This kernel is a sparse matrix, so it returns a compressed column storage
 *  form that is efficient in arithmetical operations.
 *
 *  \param min_freq Lower bound of frequency range
 *  \param max_freq Higher bound of frequency range
 *  \param sample_rate Sample rate of a source signal
 *  \param bins Number of bins per octave
 *  \param[out] height Height of a returned matrix
 *  \param[out] width Width of a returned matrix
 */
cs_ci* cq_make_kernel(int min_freq, int max_freq, int sample_rate, int bins, int* height,
                      int* width);

/**
 *  \brief Computes constant Q transform of a given data
 *
 *  \param data Source data for computation
 *  \param kernel Kernel obtained from cq_make_kernel
 *  \param height Height of a kernel
 *  \param width Width of a kernel
 *  \param plan FFTW plan for computing FFT, it has to be of R2C kind
 *  \param buffer Temporary buffer for computation, reused for optimization
 */
fftw_complex* cq_const_q_transform(double* data, const cs_ci* kernel, int height, int width,
                                   const fftw_plan plan, fftw_complex* buffer);

/**
 *  \brief Wraps constant Q transform function
 *
 *  Creates and deletes FFTW's plan and buffer for optimization purposes
 *
 *  \param data Source data for computation
 *  \param kernel Kernel obtained from cq_make_kernel
 *  \param height Height of a kernel
 *  \param width Width of a kernel
 *  \param[in] indices Array with consecutive offsets in data
 *  \param indices_size Size of indices array
 */
fftw_complex* cq_const_q_wrap(double* data, const cs_ci* kernel, int height, int width,
                              const int* indices, int indices_size);

/**
 *  \brief Computes short time constant Q transform
 *
 *  Tested only with arrays which length is a power of 2.
 *
 *  \param data Source data
 *  \param data_length Length of source data
 *  \param min_freq Lower bound of frequency range
 *  \param max_freq Higher bound of frequency range
 *  \param sample_rate Sample rate of given signal
 *  \param bins Number of bins per octave
 *  \param step Difference between consecutive slices of data used in the algorithm
 *  \param[out] height Height of a returned matrix
 *  \param[out] width Width of a returned matrix
 */
fftw_complex* cq_short_time_constq_transform(double* data, int data_length, int min_freq,
                                             int max_freq, int sample_rate, int bins, int step,
                                             int* height, int* width);

#endif /* CQ_FFT_H */