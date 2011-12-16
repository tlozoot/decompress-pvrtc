// My reimplementation of the PVRT decompression on CPU...not yet working

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "pvrtdecompress.h"

#define Y_BLOCK_SIZE 4
#define X_BLOCK_SIZE 4
#define X_BLOCK_SIZE_4BPP 4
#define X_BLOCK_SIZE_2BPP 8

typedef struct {
    uint16_t color_a;
    uint16_t color_b;
    uint32_t mod_data;
} PVRT_BLOCK_RAW;

int shiftl_32(unsigned int *x, int n);
uint8_t shiftl_16(uint16_t *x, int n);
int shiftr_32(unsigned int *x, int n);
uint8_t shiftr_16(uint16_t *x, int n);
uint8_t replicate_top_bit(uint8_t x);
void extract_color_info(RGBQUAD *pixel, PVRT_BLOCK_RAW *raw_block, int is_b);
int get_mode(PVRT_BLOCK_RAW *raw_block);
int next_weight(uint32_t *mod_data, int mode);
void print_block(RGBQUAD *rgba_block);

void pvrtdecompress(const void *input_buf, const int is_2bpp,
                    const int xmax, const int ymax,
                    unsigned char *result_buf)
{    
    int num_blocks = (xmax * ymax / 12) + 11;
    printf("%d %d\n", sizeof(PVRT_BLOCK_RAW), num_blocks);
    int input_size = sizeof(PVRT_BLOCK_RAW) * num_blocks;
    printf("decompress: %x %x\n", (unsigned) input_buf, (unsigned) input_buf + input_size);
    
    PVRT_BLOCK_RAW *input_pos = (PVRT_BLOCK_RAW *) input_buf;
    
    int offset = 0;
    
    for (int block_i = 0; block_i < num_blocks; block_i++)
    {
        PVRT_BLOCK_RAW *input_pos = (PVRT_BLOCK_RAW *) input_buf + block_i;
        printf("%x ", input_pos);
        // Get nice structs color color a and b
        RGBQUAD base_a;
        RGBQUAD base_b;
        extract_color_info(&base_a, input_pos, 0);
        extract_color_info(&base_b, input_pos, 1);
        
        int mode = get_mode(input_pos);
        
        // printf("%x ", input_pos->color_a);
        // printf("%x%x%x%x\n", base_a.a, base_a.r, base_a.g, base_a.b);
        // printf("%x%x%x%x\n", base_b.a, base_b.r, base_b.g, base_b.b);
        
        
        // Fill out the entire block of pixels based on modulation data
        RGBQUAD rgba_block[Y_BLOCK_SIZE][X_BLOCK_SIZE];
        for (int i = 0; i < Y_BLOCK_SIZE; i++) {
            for (int j = 0; j < X_BLOCK_SIZE; j++)
            {
                int weight = next_weight(&(input_pos->mod_data), mode);
                rgba_block[i][j].a = base_a.a + (weight * (base_b.a - base_a.a) >> 3);
                rgba_block[i][j].r = base_a.r + (weight * (base_b.r - base_a.r) >> 3);
                rgba_block[i][j].g = base_a.g + (weight * (base_b.g - base_a.g) >> 3);
                rgba_block[i][j].b = base_a.b + (weight * (base_b.b - base_a.b) >> 3);
                
                if ((mode == 1) && (weight == 4))
                    rgba_block[i][j].a = 0;
                
                // print_block(&rgba_block[i][j]);
            }
        }
        // printf("\n");
        
        // Write our rgba block into our result buffer
        memcpy(result_buf + offset, rgba_block, 16 * sizeof(RGBQUAD));
        // printf("%x ", result_buf + offset);
        
        input_pos += sizeof(PVRT_BLOCK_RAW);
        offset += 16 * sizeof(RGBQUAD);
    }   
    printf("\n");
}


void extract_color_info(RGBQUAD *pixel, PVRT_BLOCK_RAW *raw_block, int is_b)
{
    uint16_t raw_pixel = is_b ? raw_block->color_b : raw_block->color_a;
    int is_opaque = shiftl_16(&raw_pixel, 1);
    
    if (is_opaque) {
        pixel->r = shiftl_16(&raw_pixel, 5); 
        pixel->g = shiftl_16(&raw_pixel, 5);
        pixel->b = is_b ? replicate_top_bit(shiftl_16(&raw_pixel, 5)) : shiftl_16(&raw_pixel, 5);            
        pixel->a = 0xf;
    }
    
    else {
        pixel->a = shiftl_16(&raw_pixel, 3) << 1;
        pixel->r = replicate_top_bit(shiftl_16(&raw_pixel, 4) << 1);
        pixel->g = replicate_top_bit(shiftl_16(&raw_pixel, 4) << 1);
        pixel->b = shiftl_16(&raw_pixel, 4) << 1;
        if (is_b)
            pixel->b |= pixel->b >> 3;
        else
            pixel->b = replicate_top_bit(pixel->b);
    }
    
}

int next_weight(uint32_t *mod_data, int mode)
{
    if (mode) { // mode = 1
        switch (shiftr_32(mod_data, 2)) {        
            case 0: return 0;
            case 1: return 4;        
            case 2: return 4;
            case 3: return 8;
        }
    }
    
    else {
        switch (shiftr_32(mod_data, 2)) {        
            case 0: return 0;
            case 1: return 3;        
            case 2: return 5;
            case 3: return 8;
        }
    }
    return 0;
}

int get_mode(PVRT_BLOCK_RAW *raw_block)
{
    unsigned copy = raw_block->mod_data;
    return shiftl_32(&copy, 1);
}

// Returns the n left bits of x and shifts x by n bits
int shiftl_32(unsigned int *x, int n)
{
    int result = *x >> (32 - n);
    *x <<= n;
    return result;
}

// Returns the n left bits of x and shifts x by n bits
uint8_t shiftl_16(uint16_t *x, int n)
{
    uint8_t result = *x >> (16 - n);
    *x <<= n;
    return result;
}


// Returns the n rightmost bits of x and shifts x by n bits
int shiftr_32(unsigned int *x, int n)
{
    int result = (*x << (32 - n)) >> (32 - n);
    *x >>= n;
    return result;
}

uint8_t shiftr_16(uint16_t *x, int n)
{
    uint8_t result = (*x << (32 - n)) >> (32 - n);
    *x >>= n;
    return result;
}


void print_block(RGBQUAD *rgba_block)
{
    printf("%x", rgba_block->r);
    printf("%x", rgba_block->g);
    printf("%x", rgba_block->b);
    printf("%x", rgba_block->a);
    printf(" ");
}

uint8_t replicate_top_bit(uint8_t x)
{
    return x | x >> 4;
}

