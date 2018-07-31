/*
 * @file    main.c
 * @author  Piotr Gregor <piotr@dataandsignal.com>
 * @brief   Low pass filter.
 * @date	29 July 2018
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <arpa/inet.h>

#define BUF_SZ 400

#define SEGMENT_LEN 16

#define FILTER_LEN 257
#define FILTER_HALF_BAND 100
#define FILTER_CENTRE 0

uint16_t Fs = 8000;
double h[FILTER_LEN];
double w[FILTER_LEN];

double test_arr[16] = { 0, 0.133368129077583, 0.265840410313570, 0.396528503497067, 0.524559017724909, 0.649080821755675, 0.769272158834902, 0.884347503518051, 0.993564100303636, 1.096228126707953, 1.191700426740037, 1.279401764541372, 1.358817552205566, 1.429502010451180, 1.491081725844580, 1.543258573614740 };
double test_res[SEGMENT_LEN + FILTER_LEN - 1] = {};

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

void
usage(const char *name)
{
	fprintf(stderr, "Usage:\n");
	fprintf(stderr, "%s <input file>\n\n", name);
}

static double sinc(double x) {

	if (abs(x) < 1.0 / Fs) {
		return 1;
	}

	return sin(M_PI * x) / (M_PI * x);
}

static void hamming(uint16_t len) {

	uint16_t half = 0;
	int i = 0;

	if (len % 2) {

		// odd length window
		half = (len + 1) / 2;	// 129
	} else {
		// even length window
		half = len / 2;
	}

	while (i < half) {
		
		w[i] = 0.54 - 0.46 * cos(2 * M_PI * ((double) i / (double) (len - 1)));
		w[FILTER_LEN - 1 - i] = w[i];
		i++;
	}

}

static void create_filter(uint16_t len) {
	
	int n = 0;
	int L2 = (len - 1) / 2;
	
	
	hamming(len);
	
	while (n < len) {

		//shift_re[n] = cos(2 * PI * FILTER_CENTRE * dt * n);
		//shift_im[n] = sin(2 * PI * FILTER_CENTRE * dt * n);

		h[n] = (2.0 * FILTER_HALF_BAND / Fs) * sinc(2.0 * FILTER_HALF_BAND * (n - L2) / Fs);
		n++;
	}
	n = 0;
	while (n < len) {

		h[n] = h[n] * w[n];
		n++;
	}
}

static void conv(double *a, uint16_t a_len, double *h, uint16_t h_len, double *res) {
	
	int t = 0, tau = 0;
	int c_len = a_len + h_len - 1;
	int tau_min = 0, tau_max = 0;

	while (t < c_len) {
	
		res[t] = 0;
/**
		tau_min = (t >= h_len - 1) ? t - (h_len - 1) : 0;
		tau_max = (t < h_len - 1) ? t : h_len - 1;

		for (tau = tau_min; tau <= tau_max; tau++) {
			
			res[t] += a[tau] * h[t - tau];
		}
		**/
		
		tau_min = (t >= h_len - 1) ? t - (h_len - 1) : 0;
		tau_max = (t < h_len - 1) ? min(a_len - 1, t) : min(a_len - 1, h_len - 1);
		//tau_max = min(a_len - 1, h_len - 1);
		//tau_max = (t < h_len - 1) ? t : h_len - 1;

		//for (tau = tau_min; tau < a_len && (t - tau >= 0); tau++) {
		for (tau = tau_min; tau <= tau_max; tau++) {

			res[t] += a[tau] * h[t - tau];
		}

		t++;
	}
}

static double energy(double *a, int a_len) {

	int i = 0;
	double res = 0;

	while (i < a_len) {
		res += a[i] * a[i];
		i++;
	}

	return res;
}

static double energy_test(double *a, int a_len, double *h, int h_len) {

	double ea = 0, eaf = 0;
	double af[SEGMENT_LEN + FILTER_LEN - 1];

	conv(a, a_len, h, h_len, af);

	ea = energy(a, SEGMENT_LEN);
	eaf = energy(af, SEGMENT_LEN + FILTER_LEN - 1);

	if (ea < 0.0001) ea = 0.0001;

	fprintf(stdout, "a=%f\taf=%f,\t\ttest=%f\n", ea, eaf, eaf / ea);

	return eaf/ea;
}

void int2float(int16_t *input, double *output, int length) {
	int i = 0;
	
	for (; i < length; i++) {
		output[i] = (double) input[i] / ((double) INT16_MAX);
	}
}


int main(int argc, char *argv[])
{
	FILE    *fpin = NULL;
	char swap = '0';
	int16_t buf[BUF_SZ]= {0};
	size_t	read = 0, read_now = 0, count = 0;
	size_t i = 0;

	double a[SEGMENT_LEN] = {};


	if (argc < 2 || argc > 3) {
		fprintf(stderr, "\nProgram takes 1 argument (input file) and at most one optional (swap).\n\n");
		goto bpfhelp;
	}

	fpin = fopen(argv[1], "r");
	if (fpin == NULL) {
		fprintf(stderr, "\nError opening input file.\nPlease check the file name.\n\n");
		goto bpfhelp;
	}
	
	if (argc == 3) { 
		swap = *argv[2];
		if (swap != 's') {
			fprintf(stderr, "\nError. The last argument should be 's' if endiannes swapping required.\n\n");
			goto bpfhelp;
		}
	}

	create_filter(FILTER_LEN);
	conv(test_arr, SEGMENT_LEN, h, FILTER_LEN, test_res);

	while ((read_now = fread((void *) buf, sizeof(buf[0]), SEGMENT_LEN, fpin)) > 0) {

		/**
		 * Optionally swap endiannes.
		 */

		if (swap == 's') {

			i = 0;

			while (i < read_now) {
				buf[i] = htons(buf[i]);
				++i;
			}
		}

		read += read_now;						/* in items */

		if (read_now == SEGMENT_LEN) {

			int2float(buf, a, SEGMENT_LEN);
			energy_test(a, SEGMENT_LEN, h, FILTER_LEN);

			count += SEGMENT_LEN;
		}
	}

	fclose(fpin);

	fprintf(stderr, "\nDone. ");

	fprintf(stderr, "%zu samples read\n\n", count);

	if (swap == 's') {
		fprintf(stderr, "Endiannes swapped.\n\n");
	}


	return EXIT_SUCCESS;

bpfhelp:
	usage(argv[0]);
	return EXIT_FAILURE;
}
