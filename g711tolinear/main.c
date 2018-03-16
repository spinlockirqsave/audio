/*
 * @file    main.c
 * @author  Piotr Gregor <piotrek.gregor gmail com>
 * @brief   
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "g711.h"

#define BUF_SZ 400

int alaw2linear(int	a_val);


void
usage(const char *name)
{
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "%s <input file> <output file>\n\n", name);
}

int main(int argc, char *argv[])
{
    FILE    *fpin = NULL, *fpout = NULL;
    int16_t buf[BUF_SZ]= {0};
    int8_t val = 0;
	int16_t s = 0;
	size_t	read = 0, read_now = 0, count = 0;
	size_t items_n = 0;
	size_t i = 0;


    if (argc < 3) {
        fprintf(stderr, "Program takes 2 arguments (input and output files).\n\n");
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    fpin = fopen(argv[1], "r");
    if (fpin == NULL) {
        fprintf(stderr, "Error opening input file.\nPlease check the file name.\n\n");
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    fpout = fopen(argv[2], "w");
    if (fpout == NULL) {
        fprintf(stderr, "Error opening output file.\nPlease check the file name.\n\n");
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    while ((read_now = fread((void *) buf, sizeof(buf[0]), BUF_SZ, fpin)) > 0) {

		/* swap endiannes */
		i = 0;
/*		while (i < read_now) {
			buf[i] = ntohs(buf[i]);
			++i;
		}*/
		read += read_now;	/* in items */

		i = 0;				/* in items */

		while (i < read_now * 2) {

			val = * ((char *) buf + i);
			//s = alaw2linear(val);
			s = alaw_to_linear(val);


			if (fwrite(&s, sizeof(s), 1, fpout) != 1) {
				goto err;
			}

			++count;
			++i;
		}
    }

	fclose(fpin);
	fclose(fpout);

	fprintf(stderr, "Done: %zu samples written, %zu/%zu bytes read/written\n\n", count, read * 2, read * 4);

    return EXIT_SUCCESS;

err:

	return EXIT_FAILURE;
}
