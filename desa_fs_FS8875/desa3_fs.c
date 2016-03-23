/*
 * @file    desa3_fs.c
 * @date    20 Mar 2016
 * @brief   DESA implementations tweaked for speed.
 *          (DESA = Discrete Energy Separation Algorithm).
 *
 * @author  Piotr Gregor < piotrek.gregor gmail.com >
 *
 */
 
#include <stdio.h>
#include <stdint.h>
#include <memory.h>
#include <math.h>

#include "buffer.h"
#include "sma_buf.h"
#include "desa2_fs.h"
 
#define BLOCK       160         /* samples processed in each invocation */
#define TO_HZ(r, f) (((r) * (f)) / (2.0 * M_PI))
 
/*
 * Purpose: detect a tone using DESA-1 algorithm
 *
 * Parameters:
 *      input       pointer to input samples
 *      variance    the variance of the frequency estimates
 *
 * Return value: frequency estimate in Hz
 */
double
desa1(double *input, double *variance, uint32_t sample_rate)
{
    // detector variables
    static double diff0 = 0.0;  // delayed differences
    static double diff1 = 0.0;
    static double diff2 = 0.0;
    static double diff3 = 0.0;
    static double x1 = 0.0;     // delayed inputs
    static double x2 = 0.0;
    static double x3 = 0.0;
 
    double num; // numerator
    double den; // denominator
    double freq[BLOCK]; // frequency estimates
    int i;
 
    // calculate the frequency estimate for each sample
    for ( i = 0; i < BLOCK; i++ )
    {
        diff0 = input[i] - x1;
        num = diff2 * diff2 - diff1 * diff3
            + diff1 * diff1 - diff0 * diff2;
        den = x2 * x2 - x1 * x3;
        /* instead of
         * freq[i] = sample_rate * asin(sqrt(num/(8.0 * den))) / M_PI;
         * we store only kind of OMEGA */
        freq[i] = num/den;
        // handle errors - division by zero, square root of
        // negative number or asin of number > 1 or < -1
        if (isnan(freq[i]))
        {
            freq[i] = 0;
        } else  if (isinf(freq[i]))
        {
            freq[i] = 2000.0;
        }
 
        diff3 = diff2;
        diff2 = diff1;
        diff1 = diff0;
        x3 = x2;
        x2 = x1;
        x1 = input[i];
    }
 
    // calculate mean frequency
    double mean = 0.0;
    for ( i = 0; i < BLOCK; i++ )
    {
        mean += freq[i];
    }
    mean /= (double)BLOCK;
 
    // calculate the variance in the frequency estimates
    *variance = 0.0;
    for ( i = 0; i < BLOCK; i++ )
    {
        *variance += freq[i] * freq[i];
    }
    *variance /= (double)BLOCK;
    *variance -= (mean * mean);
 
    return mean;
}
 
/*
 * Purpose: detect a tone using DESA-2 algorithm
 *
 * Parameters:
 *      input       pointer to input samples
 *      variance    the variance of the frequency estimates
 *
 * Return value: frequency estimate in Hz
 */
double
desa2(double *input, double *variance, uint32_t sample_rate)
{
    // detector variables
    static double diff0 = 0.0;  // delayed differences
    static double diff1 = 0.0;
    static double diff2 = 0.0;
    static double x1 = 0.0;     // delayed inputs
    static double x2 = 0.0;
    static double x3 = 0.0;
 
    double num; // numerator
    double den; // denominator
    double freq[BLOCK]; // frequency estimates
    int i;
 
    // calculate the frequency estimate for each sample
    for (i = 0; i < BLOCK; i++)
    {
        diff0 = input[i] - x2;
        num = diff1 * diff1 - diff0 * diff2;
        den = x2 * x2 - x1 * x3;
        /* instead of
         * freq[i] = sample_rate * asin(sqrt(num/(4.0 * den)))
         * / (2.0 * M_PI);
         * we store only kind of OMEGA */
        freq[i] = num/den;
        // handle errors - division by zero, square root of
        // negative number or asin of number > 1 or < -1
        if (isnan(freq[i]))
        {
            freq[i] = 0;
        } else  if (isinf(freq[i]))
        {
            freq[i] = 2000.0;
        }
 
        diff2 = diff1;
        diff1 = diff0;
        x3 = x2;
        x2 = x1;
        x1 = input[i];
    }
 
    // calculate mean frequency
    double mean = 0.0;
    for (i = 0; i < BLOCK; i++)
    {
        mean += freq[i];
    }
    mean /= (double)BLOCK;
 
    // calculate the variance in the frequency estimates
    *variance = 0.0;
    for ( i = 0; i < BLOCK; i++ )
    {
        *variance += freq[i] * freq[i];
    }
    *variance /= (double)BLOCK;
    *variance -= (mean * mean);
 
    return mean;
}

/*
 * Purpose: detect a tone using DESA-2 algorithm
 *          as in FreeSWITCH's current implementation.
 *
 * Parameters:
 *      input       pointer to input samples
 *      variance    the variance of the frequency estimates
 *
 * Return value: frequency estimate in Hz
 */
/*! Number of points in desa2 sample */
#define P (5)
/*! Conversion to Hertz */
#define TO_HZ(r, f) (((r) * (f)) / (2.0 * M_PI))
void
desa2_freeswitch_int(int16_t *input, double *mean1, double *mean2, double *var1, double *var2, uint32_t sample_rate)
{
    int i;
    circ_buffer_t   b;
	sma_buffer_t    sma_b;
	sma_buffer_t    sqa_b;
    double freq[BLOCK]; // frequency estimates

    INIT_CIRC_BUFFER(&b, BLOCK);
	INIT_SMA_BUFFER(&sma_b, 160);
	INIT_SMA_BUFFER(&sqa_b, 160);

	INSERT_INT16_FRAME(&b, input, BLOCK);

    // calculate the frequency estimate for each sample as in FS
    for (i = 0; i < (BLOCK - P); i++)
    {
        freq[i] = desa2_fs_tweaked(&b, i);
        APPEND_SMA_VAL(&sma_b, freq[i]);
        APPEND_SMA_VAL(&sqa_b, freq[i] * freq[i]);
	*var1 = sqa_b.sma - (sma_b.sma * sma_b.sma);
/*	printf("<<< AVMD v[%f] f[%f][%f]Hz sma[%f][%f]Hz sqa[%f]\tsample[%d]\t[%d] >>>\n",
            *var1, freq[i], TO_HZ(sample_rate, freq[i]), sma_b.sma, TO_HZ(sample_rate, sma_b.sma), sqa_b.sma, i, input[i]); */
        printf("----Desa2_fs_int: Mean kind-of-freq = %f, var = %f, REAL_FREQ = %f\t\tsample[%f]\n",
            freq[i], *var1, TO_HZ(sample_rate, 0.5 * (double)acos(freq[i])), input[i]);
    }
    /* set mean */
    *mean1 = sma_b.sma;
    /* calculate the variance */
	*var1 = sqa_b.sma - (sma_b.sma * sma_b.sma);

    /* for comparison calculate mean2 frequency & var2 */
    double mean = 0.0;
    for (i = 0; i < (BLOCK - P); i++)
    {
        mean += freq[i];
    }
    mean /= (double) (BLOCK - P);
    *mean2 = mean;
    *var2 = 0.0;
    for (i = 0; i < (BLOCK - P); i++ )
    {
        *var2 += freq[i] * freq[i];
    }
    *var2 /= (double)(BLOCK - P);
    *var2 -= (mean * mean);
 
    free(b.buf);
    return;
}

void
desa2_freeswitch_double(double *input, double *mean1, double *mean2, double *var1, double *var2, uint32_t sample_rate)
{
    int i;
    circ_buffer_t   b;
	sma_buffer_t    sma_b;
	sma_buffer_t    sqa_b;
    double freq[BLOCK]; // frequency estimates

    INIT_CIRC_BUFFER(&b, BLOCK);
	INIT_SMA_BUFFER(&sma_b, 40);
	INIT_SMA_BUFFER(&sqa_b, 40);

	INSERT_DOUBLE_FRAME(&b, input, BLOCK);

    // calculate the frequency estimate for each sample as in FS
    for (i = 0; i < (BLOCK - P); i++)
    {
        freq[i] = desa2_fs_tweaked(&b, i);
        APPEND_SMA_VAL(&sma_b, freq[i]);
        APPEND_SMA_VAL(&sqa_b, freq[i] * freq[i]);
	*var1 = sqa_b.sma - (sma_b.sma * sma_b.sma);
/*	printf("<<< AVMD v[%f] f[%f][%f]Hz sma[%f][%f]Hz sqa[%f]\tsample[%d]\t[%f][%f]>>>\n",
            *var1, freq[i], TO_HZ(sample_rate, freq[i]), sma_b.sma, TO_HZ(sample_rate, sma_b.sma), sqa_b.sma, i, input[i], GET_SAMPLE((&b), i)); */
        printf("----Desa2_fs_double: Mean kind-of-freq = %f, var = %f, REAL_FREQ = %f\t\tsample[%f]\n",
            freq[i], *var1, TO_HZ(sample_rate, 0.5 * (double)acos(freq[i])), input[i]);
    }
    /* set mean */
    *mean1 = sma_b.sma;
    /* calculate the variance */
	*var1 = sqa_b.sma - (sma_b.sma * sma_b.sma);

    /* for comparison calculate mean2 frequency & var2 */
    double mean = 0.0;
    for (i = 0; i < (BLOCK - P); i++)
    {
        mean += freq[i];
    }
    mean /= (double) (BLOCK - P);
    *mean2 = mean;
    *var2 = 0.0;
    for (i = 0; i < (BLOCK - P); i++ )
    {
        *var2 += freq[i] * freq[i];
    }
    *var2 /= (double)(BLOCK - P);
    *var2 -= (mean * mean);
 
    free(b.buf);
    free(sma_b.data);
    free(sqa_b.data);
    return;
}

// convert 16 bit ints to doubles
void
intToFloat(int16_t *input, double *output, int length)
{
    int i;
 
    for (i = 0; i < length; i++) {
        /* samples can be divided by any values and results
         * will stay the same because our DESA estimator
         * divides differences of TKEO values, so A(mplitude)
         * is left out from the equation. What matters is
         * only the ratio of samples value */
        output[i] = (double)input[i] / ((double) 1.0);
    }
}
 
int
main(int argc, char *argv[])
{
    int16_t intData[BLOCK];
    double inputData[BLOCK];
    double frequency, freq2;
    double variance, var2;
    int numWords;
    int sampleCount, i;
    char *inFileName;
    FILE *inFile;
    uint32_t sample_rate;
 
    inFileName = NULL;
 
    if (argc == 3)
    {
        sample_rate = atoi(argv[1]);
        inFileName = argv[2];
        printf("input file name = %s\n",inFileName);
    }
    else
    {
        printf("Incorrect arguments, usage:\n\tdesa <sample_rate filename>\n");
        return(1);
    }
 
    inFile = fopen(inFileName,"rb");
    if (inFile == NULL)
    {
        printf("Exiting.Cannot open input file %s\n",inFileName);
        return(1);
    }
 
    // start counting frames
    sampleCount = 0;
 
    numWords = fread(intData, sizeof(int16_t), BLOCK, inFile );
 
    // until end of file
    while( numWords == BLOCK )
    {
        intToFloat(intData, inputData, numWords);
/*
        printf("\nframe = %d\n",sampleCount);
        for (i = 0; i < BLOCK; ++i)
        {
            printf("[%d]", intData[i]);
        }*/
       /* for(i = 0; i < 5; i++){
inputData[30 + i] = 20000.0 + i;
        }*/
        // get the frequency estimates
        frequency = desa1(inputData, &variance, sample_rate);
        printf("\nDesa1: Mean kind-of-freq = %f, var = %f, std dev = %f, REAL_FREQ = %f\n",
            frequency, variance, sqrt(variance), sample_rate * asin(sqrt(frequency/8.0)) / M_PI);
 
        frequency = desa2(inputData, &variance, sample_rate);
        printf("Desa2: Mean kind-of-freq = %f, var = %f, std dev = %f, REAL_FREQ = %f\n",
            frequency, variance, sqrt(variance), sample_rate * asin(sqrt(frequency/4.0)) / (2.0 * M_PI));

        desa2_freeswitch_int(intData, &frequency, &freq2, &variance, &var2, sample_rate);
        printf("Desa2_fs_int: Mean kind-of-freq = %f, var = %f, std dev = %f, freq2 = %f, var2 = %f, REAL_FREQ = %f\n",
            frequency, variance, sqrt(variance), freq2, var2, TO_HZ(sample_rate, 0.5 * (double)acos(frequency)));

        desa2_freeswitch_double(inputData, &frequency, &freq2, &variance, &var2, sample_rate);
        printf("Desa2_fs_double: Mean kind-of-freq = %f, var = %f, std dev = %f, freq2 = %f, var2 = %f, REAL_FREQ = %f\n",
            frequency, variance, sqrt(variance), freq2, var2, TO_HZ(sample_rate, 0.5 * (double)acos(frequency)));
 
        sampleCount += BLOCK;
        numWords = fread(intData, sizeof(int16_t), BLOCK, inFile );
    }
 
    printf("\nFinished. sampleCount = %d\n",sampleCount);
 
    fclose( inFile );
    return 0;
}
