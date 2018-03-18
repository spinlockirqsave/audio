/*
 * @file    main.c
 * @author  Piotr Gregor <piotr@dataandsignal.com>
 * @brief   Conversion from G711 A-Law/U-Law to 16 bit linear.
 * @date	18 March 2018
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "g711.h"

#define BUF_SZ 400


void
usage(const char *name)
{
	fprintf(stderr, "Usage:\n");
	fprintf(stderr, "%s <input encoding> <input file> <output file> [s]\n\n", name);
	fprintf(stderr, "For A-Law encoded input use 'a', for U-Law use 'u'.\n\nExample:\n", name);
	fprintf(stderr, "\t%s a in.alaw out.raw\n\ns is optional, put it at the end of command if endiannes should be swapped.\n\n", name);
	fprintf(stderr, "Example:\n\t%s a in.alaw out.raw s\n\n", name);
}

int main(int argc, char *argv[])
{
	FILE    *fpin = NULL, *fpout = NULL;
	char enc = 'a';
	char swap = '0';
	int16_t buf[BUF_SZ]= {0};
	int8_t val = 0;
	int16_t s = 0;
	size_t	read = 0, read_now = 0, count = 0;
	size_t items_n = 0;
	size_t i = 0;


	if (argc < 4 || argc > 5) {
		fprintf(stderr, "\nProgram takes 3 or 4 arguments (encoding, input, output file, optional swap).\n\n");
		goto g7112linhelp;
	}

	enc = *argv[1];
	if (enc != 'a' && enc != 'u') {
		fprintf(stderr, "\nError. Encoding should be 'a' for A-Law, 'u' for U-Law.\n\n");
		goto g7112linhelp;
	}

	fpin = fopen(argv[2], "r");
	if (fpin == NULL) {
		fprintf(stderr, "\nError opening input file.\nPlease check the file name.\n\n");
		goto g7112linhelp;
	}

	fpout = fopen(argv[3], "w");
	if (fpout == NULL) {
		fprintf(stderr, "\nError opening output file.\nPlease check the file name.\n\n");
		goto g7112linhelp;
	}

	if (argc == 5) { 
		swap = *argv[4];
		if (swap != 's') {
			fprintf(stderr, "\nError. The last argument should be 's' if endiannes swapping required.\n\n");
			goto g7112linhelp;
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

		while (i < read_now * 2) {

			val = * ((char *) buf + i);

			if (enc == 'a') {		
				s = alaw_to_linear(val);
			} else {
				s = ulaw_to_linear(val);
			}


			if (fwrite(&s, sizeof(s), 1, fpout) != 1) {
				goto g7112linerr;
			}

			++count;
			++i;
		}
	}

	fclose(fpin);
	fclose(fpout);

	fprintf(stderr, "\nDone. ");

	if (enc == 'a') {
		fprintf(stderr, "A-Law to 16 bit linear.\n\n");
	} else {
		fprintf(stderr, "U-Law to 16 bit linear.\n\n");
	}

	fprintf(stderr, "%zu samples written, %zu/%zu bytes read/written\n\n", count, read * 2, read * 4);

	if (swap == 's') {
		fprintf(stderr, "Endiannes swapped.\n\n");
	}


	return EXIT_SUCCESS;

g7112linhelp:
	usage(argv[0]);
	return EXIT_FAILURE;

g7112linerr:
	fprintf(stderr, "\nI/O g7112linerror while writing to %s.\n", argv[3]);
	return EXIT_FAILURE;
}
