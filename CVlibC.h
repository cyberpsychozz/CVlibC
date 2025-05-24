#pragma once 

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <ctype.h>

/*
 * CVlibC.h - Custom computer vision library in C
 * 
 * Implements basic image processing functions: convolution, median filtering,
 * Gaussian blur, Canny edge detection, grayscale conversion, etc.
 *
 * Author: Bystrykh Nikita
 * Date: May 2025
 */

// Predefined convolution kernels for various filters
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

double Emboss_kernel[9] = {
    1.0, 0.0, 0.0,
    0.0,0.0, 0.0,
    0.0, 0.0, -1.0
};

double Sharpen_kernel[9] = {
    0.0, -1.0, 0.0,
    1.0, 5.0, -1.0,
    0.0, -1.0, 0.0
};

/*
 * Compare function for quicksort.
 */
int compare_uc(const void *a, const void *b){
    unsigned char ua = *(const unsigned char*)a;
    unsigned char ub = *(const unsigned char*)b;

    return (ua > ub) - (ua < ub);
}

/*
 * Extracts a region (kernel-sized window) around a pixel at (x, y) for a specific channel.
 * Handles edge padding by filling with zeros if out-of-bounds.
 */
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

/*
 * Applies a convolution kernel to a local image area and returns the result.
 * The output is clamped between 0 and 255.
 */
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

/*
 * Applies convolution to the entire image using a specified kernel.
 * Supports optional zero padding around the edges.
 */
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

/*
 * Adds padding around the image. The padding value can be customized.
 * Used to preserve image size after convolution or filtering.
 */
unsigned char *padding(unsigned char *image, int width, int height, int channels,
    int padding_size, int filling, int *new_width_out, int *new_height_out) {
    int orig_width = width;
    int orig_height = height;

    int new_width = orig_width + 2 * padding_size;
    int new_height = orig_height + 2 * padding_size;
    *new_width_out = new_width;
    *new_height_out = new_height;

    unsigned char *result = (unsigned char *)malloc(new_width * new_height * channels);

    for (int y = 0; y < new_height; y++) {
        for (int x = 0; x < new_width; x++) {
            int idx = (y * new_width + x) * channels;
            for (int c = 0; c < channels; c++) {
            result[idx + c] = (unsigned char)filling;
            }
        }
    }

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

/*
 * Calculates the median of a kernel-sized area (used in median filtering).
 */
static int median_value(unsigned char *area, int size){
    qsort(area, size, sizeof(unsigned char), compare_uc);

    if (size % 2 == 0) {
        int value = ((int)area[size / 2 - 1] + (int)area[size / 2]) / 2;
        return value;
    } else {
        return (int)area[size / 2];
    }
}

/*
 * Applies a median filter to the image. Useful for noise reduction.
 */
unsigned char *median_filter(unsigned char *image, int width, int height, int channels, int kernel_size, int padding_size, int *new_width_out, int *new_height_out){
    
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

/*
 * Generates a 2D Gaussian kernel for smoothing and noise reduction.
 */
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

/*
 * Converts an RGB image to grayscale using luminance weights.
 */
unsigned char* rgb_to_grayscale(unsigned char* image, int width, int height, int channels, int *new_width, int *new_height, int *new_channels) {
    *new_width = width;
    *new_height = height;
    *new_channels = 1;


    int gray_size = width * height;
    unsigned char* gray_image = (unsigned char*)malloc(gray_size * sizeof(unsigned char));

    for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
        int idx = (y * width + x) * channels;
        int gray_idx = y * width + x;

        int r = image[idx];
        int g = image[idx + 1];
        int b = image[idx + 2];

        gray_image[gray_idx] = (unsigned char)(0.299 * r + 0.587 * g + 0.114 * b);
        }
    }


    return gray_image;
}

/*
 * Direction offsets used for edge tracing (8-connected neighborhood).
 */
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

/*
 * Recursively follows weak edges connected to strong edges in Canny edge detection.
 */
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

/*
 * Performs bilinear interpolation for subpixel accuracy (used in non-maximum suppression).
 */
double bilinear_interpolate(unsigned char *image, int width, int height, int channels, double x, double y, int c) {
    int x0 = (int)floor(x);
    int y0 = (int)floor(y);
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    if (x0 < 0 || x1 >= width || y0 < 0 || y1 >= height)
        return 0;

    double dx = x - x0;
    double dy = y - y0;

    double I00 = image[(y0 * width + x0) * channels + c];
    double I10 = image[y0 * width + x1];
    double I01 = image[y1 * width + x0];
    double I11 = image[y1 * width + x1];

    double I_top = I00 * (1 - dx) + I10 * dx;
    double I_bottom = I01 * (1 - dx) + I11 * dx;

    return I_top * (1 - dy) + I_bottom * dy;
}

/*
 * Implements the full Canny edge detection pipeline:
 * - Gaussian blur
 * - Grayscale conversion
 * - Gradient computation (Sobel)
 * - Non-maximum suppression
 * - Double thresholding
 * - Edge tracking by hysteresis
 */
unsigned char *Canny_Edge_detector(unsigned char *image, int width, int height, int channels, int paddig_size, int *new_width_out, int *new_height_out, int *new_channels, double upper_threshold, double lower_threshold){    

    //Step 1: Applying 5x5 Gaussian kernel to the image to remove the noise
    int kernel_size = 5;
    double sigma = 3.0;
    double *gaussian_kernel = generate_gaussian_kernel(kernel_size,sigma);
    unsigned char *blurred = convolution(image, width, height, channels, gaussian_kernel, kernel_size, paddig_size, new_width_out, new_height_out);
    int new_width = *new_width_out; 
    int new_height = *new_height_out;

    // Step 2: Turning image from RGB to Gray
    int a,b;
    unsigned char *gray = rgb_to_grayscale(blurred, width, height, channels, &a, &b, new_channels);
    free(blurred);

    // Step 3: Finding intensity gradients(filtering with Sobel kernel)
    int neww;
    int newh;
    unsigned char *Gx = convolution(gray, new_width, new_height, 1, Sobel_kernel_X, 3, 1, &new_width, &new_height);
    unsigned char *Gy = convolution (gray, new_width, new_height, 1, Sobel_kernel_Y, 3, 1, &new_width, &new_height);
    
    unsigned char *G = (unsigned char*)malloc(new_width * new_height * sizeof(unsigned char));
    double *direction = (double*)malloc(new_width * new_height * sizeof(double));
    
    for(int i = 0; i < new_width * new_height; i++){
        double gx = (double)Gx[i];
        double gy = (double)Gy[i];
        G[i] = (unsigned char)fmin(255, hypot(gx, gy));        
        direction[i] = atan2(gy, gx);
    }
    free(Gx);
    free(Gy);

    // Step 4: Nonmaximum Suppresion
    unsigned char *supressed = (unsigned char*)malloc(new_height * new_width * sizeof(unsigned char));

    for(int y = 0; y < new_height -1; y++){
        for(int x = 0; x < new_width -1; x++){

            int idx = y * new_width + x;
            double angle = direction[idx];
            unsigned char curr = G[idx];

            angle = fmod(angle + M_PI, M_PI);
            angle = angle * 180 / M_PI;

            double dx = cos(angle);
            double dy = sin(angle);

            double x1 = x + dx;
            double y1 = y + dy;
            double x2 = x - dx;
            double y2 = y - dy;

            double q = bilinear_interpolate (G, new_width, new_height, *new_channels, x1, y1, 0);
            double r = bilinear_interpolate (G, new_width, new_height, *new_channels, x2, y2, 0);

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

/*
 * Computes the new image dimensions after an affine transformation.
 */
void compute_new_size_affine(int width, int height,
                             double a, double b, double c, double d, double tx, double ty,
                             int *new_width, int *new_height,
                             double *min_x_out, double *min_y_out) {
    double corners[4][2] = {
        {0, 0},
        {width, 0},
        {0, height},
        {width, height}
    };

    double min_x = 1e9, max_x = -1e9;
    double min_y = 1e9, max_y = -1e9;

    for (int i = 0; i < 4; i++) {
        double x = corners[i][0];
        double y = corners[i][1];

        double tx_ = a * x + b * y + tx;
        double ty_ = c * x + d * y + ty;

        if (tx_ < min_x) min_x = tx_;
        if (tx_ > max_x) max_x = tx_;
        if (ty_ < min_y) min_y = ty_;
        if (ty_ > max_y) max_y = ty_;
    }

    *new_width = (int)ceil(max_x - min_x);
    *new_height = (int)ceil(max_y - min_y);
    *min_x_out = min_x;
    *min_y_out = min_y;
}

// Affine transformation matrix:
// | a  b  tx |
// | c  d  ty |

/*
 * Applies an affine transformation (e.g., rotation, scaling, translation) to the image.
 * Uses bilinear interpolation for pixel value computation.
 */
unsigned char *affine_transform(unsigned char *img, int width, int height, int channels,
                                double a, double b, double c, double d, double tx, double ty,
                                int *new_width, int *new_height, int type) {

    double min_x, min_y;
    compute_new_size_affine(width, height, a, b, c, d, tx, ty, new_width, new_height, &min_x, &min_y);

    unsigned char *result = (unsigned char*)malloc((*new_width) * (*new_height) * channels);
    if (!result) return NULL;

    for (int y = 0; y < *new_height; y++) {
        for (int x = 0; x < *new_width; x++) {

            double src_x = a * x + b * y + tx - min_x;
            double src_y = c * x + d * y + ty - min_y;

            for (int ch = 0; ch < channels; ch++) {
                double val = bilinear_interpolate(img, width, height, channels, src_x, src_y, ch);
                if (val < 0) val = 0;
                if (val > 255) val = 255;

                result[(y * (*new_width) + x) * channels + ch] = (unsigned char)val;
            }
        }
    }

    return result;
}