#pragma once 

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 


/*TODO
    Convolution - completed
    Padding
    Median filter
	Gaussian filter
	Edge detection
    smth else from Opencv
    Image resizing
    Binarisation
*/

// Функция для перемножения области изображения и ядра
//Function for multyplying area of image with kernel
static int conv_matrix_multiply(unsigned char *area, int *kernel, int size) {
    int sum = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            sum += area[i * size + j] * kernel[i * size + j];
        }
    }

    if (sum < 0) sum = 0;
    if (sum > 255) sum = 255;
    return sum;
}

// Функция для выделения области вокруг пикселя (x, y) для одного канала
static void create_area(unsigned char *image, int width, int height, int channels,
    int x, int y, int channel, int kernel_size, unsigned char *area) {
    
    int half = kernel_size / 2;

    for (int ky = 0; ky < kernel_size; ky++) {
        for (int kx = 0; kx < kernel_size; kx++) {

            int img_x = x + (kx - half);
            int img_y = y + (ky - half);

            if (img_x >= 0 && img_x < width && img_y >= 0 && img_y < height) {
                int idx = (img_y * width + img_x) * channels + channel;
                area[ky * kernel_size + kx] = image[idx];
            } else {
                area[ky * kernel_size + kx] = 0;
            }
        }
    }
}


unsigned char *convolution(unsigned char *image, int width, int height, int channels, 
                int *kernel, int kernel_size, int padding, int *new_width_out, int *new_height_out) {
    int new_width = (width - kernel_size + 2 * padding) + 1;
    int new_height = (height - kernel_size + 2 * padding) + 1;
    *new_width_out = new_width;
    *new_height_out = new_height;
    // int beg_x = padding;
    // int beg_y = padding;
    // int end_x = beg_x + (width - kernel_size + 1);
    // int end_y = beg_y + (height - kernel_size + 1);

    unsigned char *result = (unsigned char*)malloc(new_height * new_width * channels * sizeof(unsigned char));

    unsigned char *area = (unsigned char*)malloc(kernel_size * kernel_size * sizeof(unsigned char));

    
    int half = kernel_size / 2;
    for (int y = 0; y < new_height; y++) {
        for (int x = 0; x < new_width; x++) {
            
            int src_x = x + half - padding;
            int src_y = y + half - padding;

            int result_idx = (y * new_width + x) * channels;
            
            for (int c = 0; c < channels; c++) {
                create_area(image, width, height, channels, src_x, src_y, c, kernel_size, area);
                result[result_idx + c] = (unsigned char)conv_matrix_multiply(area, kernel, kernel_size);
            }
        }
    }

    free(area);
    return result;
}

unsigned char *padding (unsigned char *image, int width, int height, int channels, 
    int padding_size, int filling, int *new_width_out, int *new_height_out){
    
        width += 2 * padding_size;
    height += 2 * padding_size;
    int new_width = width;
    int new_height = height;
    *new_width_out = new_width;
    *new_height_out = new_height;

    unsigned char *result = (unsigned char*)malloc(height * width * channels * sizeof(unsigned char));
   
    // int beg_img = padding_size

    for (int y = 0; y < new_height; y++) {
        for (int x = 0; x < new_width; x++) {
            int idx = (y * new_width + x) * channels;
            for (int c = 0; c < channels; c++) {
                result[idx + c] = (unsigned char)filling;
            }
        }
    }

    
    for (int y = padding_size; y < (new_height - padding_size); y++) {
        for (int x = padding_size; x < (new_width - padding_size); x++) {
            int orig_x = x - padding_size;
            int orig_y = y - padding_size;
            int orig_idx = (orig_y * width + orig_x) * channels;
            int idx = (y * new_width + x) * channels;
            for (int c = 0; c < channels; c++) {
                result[idx + c] = image[orig_idx + c];
            }
        }
    }

    return result;
}