/*
 * @file    main.c
 * @author  Piotr Gregor <piotrek.gregor gmail com>
 * @brief   Blob detection.
*/


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


size_t
msb_pos(uint32_t u)
{
    size_t bits = 0;
    while (u != 0)
    {
        ++bits;
        u >>= 1;
    };
    return bits;
}

unsigned char *pImage;
uint32_t width = 10;
uint32_t length = 10;
uint32_t yellow = 145;
uint32_t red = 234;

int
main(void)
{
    int i;
    uint32_t size;
    unsigned char areas[256] = { 0 };

    if (msb_pos(width) + msb_pos(length) > 32)
    {
        perror("Width and/or length too big, "
                "would overflow\n");
        return -1;
    }
    size = width * length;
    pImage = malloc(size);
    if (pImage == NULL)
        return -1;
    i = 0;
    for (; i < size; ++i)
    {
        areas[pImage[i]]++;
    }
    free(pImage);
    return 0;
}
