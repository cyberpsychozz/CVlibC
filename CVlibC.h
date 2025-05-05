#pragma once 

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 

double Sobel_kernel_Y[9] = {
    1.0, 2.0, 1.0,
    0.0, 0.0, 0.0,
    -1.0, -2.0, -1.0
};

double Sobel_kernel_X[9] = {
    -1.0, 0.0, 1.0,
    -2.0, 0.0, 2.0,
    -1.0, 0.0, 1.0
};

/*TODO
    Convolution - completed
    Padding - Completed
    Median filter - Completed
	Gaussian filter(Function, that creates Gauss Kernel) - Completed
    Canny Edge detection - completed
    RGB to Gray - completed


    smth else from Opencv
    Image resizing
    Binarisation
    Threshold parameters in Canny filter
    User friendly variant of main
*/



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

//Function for multyplying area of image with kernel
static int conv_matrix_multiply(unsigned char *area, double *kernel, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            sum += area[i * size + j] * kernel[i * size + j];
        }
    }

    if (sum < 0.0) sum = 0;
    if (sum > 255.0) sum = 255;
    return sum;
}

unsigned char *convolution(unsigned char *image, int width, int height, int channels, 
    double *kernel, int kernel_size, int padding_size, int *new_width_out, int *new_height_out) {
        int new_width = (width - kernel_size + 2 * padding_size) + 1;
        int new_height = (height - kernel_size + 2 * padding_size) + 1;
        *new_width_out = new_width;
        *new_height_out = new_height;

        unsigned char *result = (unsigned char*)malloc(new_height * new_width * channels * sizeof(unsigned char));

        unsigned char *area = (unsigned char*)malloc(kernel_size * kernel_size * sizeof(unsigned char));


        int half = kernel_size / 2;
        for (int y = 0; y < new_height; y++) {
            for (int x = 0; x < new_width; x++) {

                int src_x = x + half - padding_size;
                int src_y = y + half - padding_size;

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

unsigned char *padding(unsigned char *image, int width, int height, int channels,
    int padding_size, int filling, int *new_width_out, int *new_height_out) {
    int orig_width = width;
    int orig_height = height;

    int new_width = orig_width + 2 * padding_size;
    int new_height = orig_height + 2 * padding_size;
    *new_width_out = new_width;
    *new_height_out = new_height;

    unsigned char *result = (unsigned char *)malloc(new_width * new_height * channels);

    // 
    for (int y = 0; y < new_height; y++) {
        for (int x = 0; x < new_width; x++) {
            int idx = (y * new_width + x) * channels;
            for (int c = 0; c < channels; c++) {
            result[idx + c] = (unsigned char)filling;
            }
        }
    }

    // image copying
    for (int y = 0; y < orig_height; y++) {
        for (int x = 0; x < orig_width; x++) {
            int orig_idx = (y * orig_width + x) * channels;
            int idx = ((y + padding_size) * new_width + (x + padding_size)) * channels;
            for (int c = 0; c < channels; c++) {
            result[idx + c] = image[orig_idx + c];
            }
        }
    }

    return result;
}

int compare_uc(const void *a, const void *b){
    unsigned char ua = *(const unsigned char*)a;
    unsigned char ub = *(const unsigned char*)b;

    return (ua > ub) - (ua < ub);
}

static int median_value(unsigned char *area, int size){
    qsort(area, size, sizeof(unsigned char), compare_uc);

    if (size % 2 == 0) {
        int value = ((int)area[size / 2 - 1] + (int)area[size / 2]) / 2;
        return value;
    } else {
        return (int)area[size / 2];
    }
}

unsigned char *median_filter(unsigned char *image, int width, int height, int channels, int kernel_size, int padding_size, int *new_width_out, int *new_height_out){
    
    //New parameters of image width and height
    int new_width = (width - kernel_size + 2 * padding_size) + 1;
    int new_height = (height - kernel_size + 2 * padding_size) + 1;
    *new_width_out = new_width;
    *new_height_out = new_height;

    unsigned char *result = (unsigned char*)malloc(new_height * new_width * channels * sizeof(unsigned char));

    unsigned char *area = (unsigned char*)malloc(kernel_size * kernel_size * sizeof(unsigned char));

    int half = kernel_size / 2;
    for (int y = 0; y < new_height; y++) {
        for (int x = 0; x < new_width; x++) {
            
            int src_x = x + half - padding_size;
            int src_y = y + half - padding_size;

            int result_idx = (y * new_width + x) * channels;
            
            for (int c = 0; c < channels; c++) {
                create_area(image, width, height, channels, src_x, src_y, c, kernel_size, area);
                result[result_idx + c] = (unsigned char)median_value(area, kernel_size * kernel_size);

            }
        }
    }

    free(area);
    return result;
}

double *generate_gaussian_kernel(int size, double sigma) {
    
    double *kernel = (double*)malloc(size * size * sizeof(double));

    int half = size / 2;
    double sum = 0.0;

    for (int y = -half; y <= half; y++) {
        for (int x = -half; x <= half; x++) {
            double exponent = -(x * x + y * y) / (2 * sigma * sigma);
            double value = exp(exponent) / (2 * M_PI * sigma * sigma);
            kernel[(y + half) * size + (x + half)] = value;
            sum += value;
        }
    }

    
    for (int i = 0; i < size * size; i++) {
        kernel[i] /= sum;
    }

    return kernel;
}

// Function turning 
unsigned char* rgb_to_grayscale(unsigned char* image, int width, int height, int channels, int *new_width, int *new_height, int *new_channels) {
    *new_width = width;
    *new_height = height;
    *new_channels = 1;


    int gray_size = width * height;
    unsigned char* gray_image = (unsigned char*)malloc(gray_size * sizeof(unsigned char));

    for (int i = 0; i < width * height; i++) {
        int r = image[i * channels];
        int g = image[i * channels + 1];
        int b = image[i * channels + 2];

        // Y = 0.299 * R + 0.587 * G + 0.114 * B
        gray_image[i] = (unsigned char)(0.299 * r + 0.587 * g + 0.114 * b);
    }

    return gray_image;
}

int d[16] = {
    1,  0,
    1,  1,
    0,  1,
   -1,  1,
   -1,  0,
   -1, -1,
    0, -1,
    1, -1
};

//Function for recursive following edge

void edge_following(unsigned char *input, unsigned char *output, int width, int height, int x, int y){
   for(int i = 0; i < 16; i += 2){
       int nx = x + d[i];
       int ny = y + d[i + 1];

       if (nx < 0 || nx >= width || ny < 0 || ny >= height) continue;

       int idx = ny * width + nx;

       if (input[idx] == 100 && output[idx] != 255) {
           output[idx] = 255;
           edge_following(input, output, width, height, nx, ny);
       }
   }
}



unsigned char *Canny_Edge_detector(unsigned char *image, int width, int height, int channels, int paddig_size, int *new_width_out, int *new_height_out, int *new_channels, double upper_threshold, double lower_threshold){

    // Step 1: Turning image from RGB to Gray
    int a,b;
    unsigned char *gray = rgb_to_grayscale(image, width, height, channels, &a, &b, new_channels);
    

    //Step 2: Applying 5x5 Gaussian kernel to the image to remove the noise
    int kernel_size = 5;
    double sigma = 2.0;
    double *gaussian_kernel = generate_gaussian_kernel(kernel_size,sigma);
    unsigned char *blurred = convolution(gray, width, height, 1, gaussian_kernel, kernel_size, paddig_size, new_width_out, new_height_out);
    int new_width = *new_width_out; 
    int new_height = *new_height_out;

    // Step 3: Finding intensity gradients(filtering with Sobel kernel)
    int neww;
    int newh;
    unsigned char *Gx = convolution(blurred, new_width, new_height, 1, Sobel_kernel_X, 3, 1, &new_width, &new_height);
    unsigned char *Gy = convolution (blurred, new_width, new_height, 1, Sobel_kernel_Y, 3, 1, &new_width, &new_height);
    free(blurred);

    unsigned char *G = (unsigned char*)malloc(new_width * new_height * sizeof(unsigned char));
    double *direction = (double*)malloc(new_width * new_height * sizeof(double));
    
    for(int i = 0; i < new_width * new_height; i++){
        double gx = (double)Gx[i];
        double gy = (double)Gy[i];
        G[i] = (unsigned char)fmin(255, hypot(gx, gy));        
        direction[i] = atan2(gx, gy);
    }
    free(Gx);
    free(Gy);


    // Step 4: Nonmaximum Suppresion

    unsigned char *supressed = (unsigned char*)malloc(new_height * new_width * sizeof(unsigned char));


    for(int y = 0; y < new_height - 1; y++){
        for(int x = 0; x < new_width - 1; x++){

            int idx = y * new_width + x;
            double angle = direction[idx];

            unsigned char curr = G[idx];
            unsigned char q = 0, r = 0;


            angle = fmod(angle + M_PI, M_PI);
            angle = angle * 180 / M_PI;

            if ((angle >= 0 && angle < 22.5) || (angle >= 157.5 && angle <= 180)){
                q = G[y * new_width + (x - 1)];
                r = G[y * new_width + (x + 1)];
            } else if(angle >= 22.5 && angle < 67.5){
                q = G[(y + 1) * new_width + (x - 1)];
                r = G[(y - 1) * new_width + (x + 1)];
            } else if(angle >= 67.5 && angle < 112.5){
                q = G[(y - 1) * new_width + x];
                r = G[(y + 1) * new_width + x];
            } else if(angle >= 112.5 && angle < 157.5){
                q = G[(y - 1) * new_width + (x - 1)];
                r = G[(y + 1) * new_width + (x + 1)];
            }

            if(q <=curr && r <= curr){
                supressed[idx] = curr;
            }else {
                supressed[idx] = 0;
            }
        }
    }

    free(G);
    free(direction);

    //Step 5: Double threshold 

    for(int i = 0; i < new_width * new_height; i++){
        if(supressed[i] < lower_threshold){
            supressed[i] = 0;
        }else if(supressed[i] >= lower_threshold && supressed[i] < upper_threshold){
            supressed[i] = 100;
        }else if(supressed[i] >= upper_threshold){
            supressed[i] = 255;
        }
    }

    //Step 6:Edge tracing by hysteresis
    unsigned char *result = (unsigned char*)calloc(new_width * new_height, sizeof(unsigned char));

    // int d{8} = [0, 1, -1];

    for(int y = 0; y < new_height; y++){
        for(int x = 0; x < new_width; x++){
            if (supressed[y * new_width + x] == 255){
                result[y * new_width + x] = 255;
                edge_following(supressed,result, new_width, new_height, x, y);
            }
        }
    }

    free(supressed);
    return result;
}