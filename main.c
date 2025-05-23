#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "CVlibC.h"


/*TODO
    Convolution 
    Median filter
	Gaussian filter
	Edge detection
    smth else from Opencv
    Image resizing
    Binarisation
*/

/*
    main.c:
    including header with functions 
    main gets argv from console and calling functions 
*/



int main(int argc, char *argv[]) {
    
    if(argc < 4){
        fprintf(stderr, "Not enough parameters in input!");
    }
   
    //Parameters from input
    char *operation;
    char *path_input;
    char *path_output;
    char *additional_parameters;

    printf("%d\n", argc);

    if(argc == 4){
        path_input = argv[1];
        operation = argv[2];
        path_output = argv[3];
    }else if(argc == 5){
        path_input = argv[1];
        operation = argv[2];
        additional_parameters = argv[3];
        path_output = argv[4];
    }

    printf("%s\n", path_input);
    printf("%s\n", path_output);


    int height, width, channels;
    unsigned char *img = stbi_load(path_input, &width, &height, &channels, 0);
    if (img) {
        printf("Image loaded: %dx%d, %d channels\n", width, height, channels);
    }

    double kernel[9] = {
        1.0, 0.0, 0.0,
        0.0,0.0, 0.0,
        0.0, 0.0, -1.0
    };

    int kernel_size = 3;
    int padding_size = 0;
    int filling = 128;

    int new_width, new_height, new_channels;

    unsigned char *result; 

    if(strcmp(operation, "-conv") == 0){
        printf("conv\n");
        result = convolution(img, width, height, channels, kernel, kernel_size, padding_size, &new_width, &new_height);
        new_channels = channels;
    }else if(strcmp(operation, "-padding") == 0){
        int pad_size = 8;
        printf("padding\n");
        result = padding(img, width, height, channels, pad_size, filling, &new_width, &new_height);
        new_channels = channels;
    }else if(strcmp(operation, "-median")==0){
        printf("median\n");
        result = median_filter(img, width, height, channels, kernel_size, padding_size, &new_width, &new_height);
        new_channels = channels;
    }else if(strcmp(operation, "-gauss")==0){
        double sigma = 1.0;
        double *kernel1 = (double*)generate_gaussian_kernel(kernel_size, sigma);
        printf("gauss\n");
        result = convolution(img, width, height, channels, kernel1, kernel_size, padding_size, &new_width, &new_height);
        free(kernel1);
        new_channels = channels;
    }else if(strcmp(operation, "-canny") == 0){
        
        result = Canny_Edge_detector(img, width, height, channels, padding_size, &new_width, &new_height, &new_channels, up_thres, low_thres);
        // result = rgb_to_grayscale(img, width, height, channels, &new_width, &new_height, &new_channels);
        // channels = 1;
        // new_width = width;
        // new_height = height;
    }

    // stbi_write_png(path_output, new_width, new_height, channels, result, 100);
    stbi_write_jpg(path_output,new_width, new_height, new_channels, result, 100);

    if(result){
        printf("image writen\n");
    }
    
    // printf("%ld\n",strlen(img));
    stbi_image_free(img);
    free(result);
    return 0;
}