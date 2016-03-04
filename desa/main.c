/*
 * @file    main.c
 * @date    26 Feb 2016
 * @brief   Tone detection using the Teager Kaiser Energy Operator to
 *          measure frequency (Desa-1 and Desa-2 algorithms)
 *          (DESA = Discrete Energy Separation Algorithm).
 *
 * Copyright 2011 Shawn Stevenson
 *
 * Modifications:   Piotr Gregor <piotrek . gregor at gmail.com>
 *
 */
 
#include <stdio.h>
#include <stdint.h>
#include <memory.h>
#include <math.h>
 
#define BLOCK   160          // samples processed each invocation
#define SAMPLE_RATE 8000.0  // sample rate of input signal
 
/*
 * Purpose: detect a tone using DESA-1 algorithm
 *
 * Parameters:
 *      input       pointer to input samples
 *      variance    the variance of the frequency estimates
 *
 * Return value: frequency estimate in Hz
 */
double desa1( double *input, double *variance )
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
        freq[i] = SAMPLE_RATE * asin(sqrt(num/(8.0 * den))) / M_PI;
        // handle errors - division by zero, square root of
        // negative number or asin of number > 1 or < -1
        if (isnan(freq[i]))
        {
            freq[i] = 0;
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
double desa2( double *input, double *variance )
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
    for ( i = 0; i < BLOCK; i++ )
    {
        diff0 = input[i] - x2;
        num = diff1 * diff1 - diff0 * diff2;
        den = x2 * x2 - x1 * x3;
        freq[i] = SAMPLE_RATE * asin(sqrt(num/(4.0 * den)))
            / (2.0 * M_PI);
        // handle errors - division by zero, square root of
        // negative number or asin of number > 1 or < -1
        if (isnan(freq[i]))
        {
            freq[i] = 0;
        }
 
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
 
// convert 16 bit ints to doubles
void intToFloat( int16_t *input, double *output, int length )
{
    int i;
 
    for ( i = 0; i < length; i++ ) {
        output[i] = (double)input[i];
    }
}
 
int main( int argc, char *argv[] )
{
    int16_t intData[BLOCK];
    double inputData[BLOCK];
    double frequency;
    double variance;
    int numWords;
    int sampleCount, i;
    char *inFileName;
    FILE *inFile;
 
    inFileName = NULL;
 
    if ( argc == 2 )
    {
        inFileName = argv[1];
        printf("input file name = %s\n",inFileName);
    }
    else
    {
        printf("Incorrect arguments, usage: desa filename\n");
        return(1);
    }
 
    inFile = fopen(inFileName,"rb");
    if ( !inFile )
    {
        printf("Exiting.Cannot open input file %s\n",inFileName);
        fclose( inFile );
        return(1);
    }
 
    // start counting frames
    sampleCount = 0;
 
    numWords = fread(intData, sizeof(int16_t), BLOCK, inFile );
 
    // until end of file
    while( numWords == BLOCK )
    {
        intToFloat( intData, inputData, numWords );
        printf("\nsample = %d\n",sampleCount);
        for (i = 0; i < BLOCK; ++i)
        {
            printf("[%d]", intData[i]);
        }
 
        // get the frequency estimates
        frequency = desa1( inputData, &variance );
        printf("\nDesa1: Mean freq = %f, var = %f, std dev = %f",
            frequency, variance, sqrt(variance));
 
        frequency = desa2( inputData, &variance );
        printf("\nDesa2: Mean freq = %f, var = %f, std dev = %f\n",
            frequency, variance, sqrt(variance));
 
        sampleCount += BLOCK;
        numWords = fread(intData, sizeof(int16_t), BLOCK, inFile );
    }
 
    printf("\nFinished. sampleCount = %d\n",sampleCount);
 
    fclose( inFile );
    return 0;
}
