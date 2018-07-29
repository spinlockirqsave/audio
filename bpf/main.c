/*
 * @file    main.c
 * @author  Piotr Gregor <piotr@dataandsignal.com>
 * @brief   Bandpass filter.
 * @date	29 July 2018
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <float.h>

#define BUF_SZ 400


void
usage(const char *name)
{
	fprintf(stderr, "Usage:\n");
	fprintf(stderr, "%s <input file>\n\n", name);
}

int main(int argc, char *argv[])
{
	FILE    *fpin = NULL, *fpout = NULL;
	char enc = 'a';
	char swap = '0';
	int16_t buf[BUF_SZ]= {0};
	int16_t val = 0;
	double s = 0;
	size_t	read = 0, read_now = 0, count = 0;
	size_t items_n = 0;
	size_t i = 0;


	if (argc < 3 || argc > 4) {
		fprintf(stderr, "\nProgram takes 2 arguments (input and output file) and at most one optional (swap).\n\n");
		goto bpfhelp;
	}

	fpin = fopen(argv[1], "r");
	if (fpin == NULL) {
		fprintf(stderr, "\nError opening input file.\nPlease check the file name.\n\n");
		goto bpfhelp;
	}

	fpout = fopen(argv[2], "w");
	if (fpout == NULL) {
		fprintf(stderr, "\nError opening output file.\nPlease check the file name.\n\n");
		goto bpfhelp;
	}
	
	if (argc == 4) { 
		swap = *argv[3];
		if (swap != 's') {
			fprintf(stderr, "\nError. The last argument should be 's' if endiannes swapping required.\n\n");
			goto bpfhelp;
		}
	}

	while ((read_now = fread((void *) buf, sizeof(buf[0]), BUF_SZ, fpin)) > 0) {

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

		i = 0;									/* in items */

		while (i < read_now) {

			val = buf[i];

			s = (double) val / INT16_MAX;
	
			fprintf(stderr, "Sample read/written: %d\t%f\n", val, s);


			if (fwrite(&s, sizeof(s), 1, fpout) != 1) {
				goto bpferr;
			}

			++count;
			++i;
		}
	}

	fclose(fpin);
	fclose(fpout);

	fprintf(stderr, "\nDone. ");

	fprintf(stderr, "%zu samples written, %zu bytes read/written\n\n", count, read * 2);

	if (swap == 's') {
		fprintf(stderr, "Endiannes swapped.\n\n");
	}


	return EXIT_SUCCESS;

bpfhelp:
	usage(argv[0]);
	return EXIT_FAILURE;

bpferr:
	fprintf(stderr, "\nI/O eror while writing to %s.\n", argv[2]);
	return EXIT_FAILURE;
}
