#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "CVlibC.h"

void print_menu() {
    printf("\nChoose operation:\n");
    printf("1 - Convolution\n");
    printf("2 - Padding \n");
    printf("3 - Median filter\n");
    printf("4 - Gaussian filter\n");
    printf("5 - Image embossing\n");
    printf("6 - Image sharpenning\n");
    printf("7 - Converting an image to grayscale\n");
    printf("8 - Canny edge detector\n");
    printf("9 - Affine transformation\n");
    printf("Your choice: ");
}



char *get_file_extension(const char *filename) {
    char *dot = strrchr(filename, '.');
    if (!dot || dot == filename) return "";
    return dot + 1;
}

int main() {
    char input_path[256];
    char output_path[256];
    int operation_choice;
    int kernel_size = 3;
    int padding_size = 0;
    int fill_value = 255;
    double sigma = 1.0;
    int affine_transform_type = 0;

    
    double kernel[9] = {
        1.0, 0.0, 0.0,
        0.0,0.0, 0.0,
        0.0, 0.0, -1.0
    };

    // 1. Getting the input file name
    printf("Enter path to input image file: ");
    if (scanf("%255s", input_path) != 1) {
        fprintf(stderr, "Invalid input.\n");
        return 1;
    }

    

    // 2. Image loading
    int width, height, channels;
    unsigned char *img = stbi_load(input_path, &width, &height, &channels, 0);
    if (!img) {
        fprintf(stderr, "Image loading failed.\n");
        return 1;
    }
    printf("Image loaded: %dx%d, %d channels\n", width, height, channels);
    int new_width, new_height, new_channels = channels;

    // 3. Operation choice menu
    print_menu();
    if (scanf("%d", &operation_choice) != 1) {
        fprintf(stderr, "Invalid operation input.\n");
        return 1;
    }


    

    //4. Optional image padding
    printf("Do You want to pad image? Y/N\n");
    char choice1[1];
    if (scanf("%1s", choice1) != 1) {
        fprintf(stderr, "Invalid input.\n");
        return 1;
    }


    if(strcmp(choice1, "Y") == 0){
        printf("Enter padding size:\n");
        if(scanf("%d", &padding_size) != 1){
            fprintf(stderr, "Invalid input.\n");
            return 1;
        }
        img = padding(img, width, height, channels, padding_size, fill_value, &new_width, &new_height);
    }



    // 5. Adittional parameters
    if (operation_choice == 1) {
        printf("Do You want to use your own kernel? Y/N\n");
        char choice[1];
        if (scanf("%1s", choice) != 1) {
            fprintf(stderr, "Invalid input.\n");
            return 1;
        }
        if (strcmp(choice, "Y") == 0){
            printf("Enter kernel size (3 or 5 for example): ");
            if (scanf("%d", &kernel_size) != 1) {
                fprintf(stderr, "Invalid input.\n");
                return 1;
            }

            printf("Enter kernel:\n");
            for(int i = 0; i < kernel_size; i++){
                for(int j = 0; j < kernel_size; j++){
                    if (scanf("%lf", &kernel[i * kernel_size + j]) != 1) {
                        fprintf(stderr, "Invalid kernel value.\n");
                        return 1;
                    }
                }
            }
        }

    }
    if (operation_choice != 7 && !(strcmp(choice1, "Y") == 0)) {
    printf("Enter padding size (1 for example): ");
    if (scanf("%d", &padding_size) != 1) {
        fprintf(stderr, "Invalid padding.\n");
        return 1;
    }
}
if (operation_choice == 2) {
    printf("Enter filling value (0 for example): ");
    if (scanf("%d", &fill_value) != 1) {
        fprintf(stderr, "Invalid fill value.\n");
        return 1;
    }
}
if (operation_choice == 4) {
    printf("Enter sigma value for Gaussian filter (1.0 for example): ");
    if (scanf("%lf", &sigma) != 1) {
        fprintf(stderr, "Invalid sigma.\n");
        return 1;
    }
}
if (operation_choice == 9){
    printf("What type of Affine transformation do you want:\n 1 - Scaling\n 2 - Rotating\n 3 - Translation\n 4 - Shearing\n");
    if (scanf("%d", &affine_transform_type) != 1 || affine_transform_type < 1 || affine_transform_type > 4) {
        fprintf(stderr, "Invalid affine transformation type.\n");
        return 1;
    }
}


    // 6. Image processing
    unsigned char *result = NULL;
    

    switch (operation_choice) {
        case 1: {
            result = convolution(img, width, height, channels, kernel, kernel_size, padding_size, &new_width, &new_height);
            break;
        }
        case 2:
            result = padding(img, width, height, channels, padding_size, fill_value, &new_width, &new_height);
            break;
        case 3:
            result = median_filter(img, width, height, channels, kernel_size, padding_size, &new_width, &new_height);
            break;
        case 4: {
            double *gauss_kernel = generate_gaussian_kernel(kernel_size, sigma);
            result = convolution(img, width, height, channels, gauss_kernel, kernel_size, padding_size, &new_width, &new_height);
            free(gauss_kernel);
            break;
        }
        case 5:
            result = convolution(img, width, height, channels, Emboss_kernel, 3, padding_size, &new_width, &new_height);
            break;
        case 6:{
            result = convolution(img, width, height, channels, Sharpen_kernel, 3, padding_size, &new_width, &new_height);
            break;
        }
        case 7:{
            result = rgb_to_grayscale(img, width, height, channels, &new_width, &new_height, &new_channels);
            break;
        }
        case 8:{
            double up_thres = 255 * 0.25;
            double low_thres = 255 * 0.15;
            result = Canny_Edge_detector(img, width, height, channels, padding_size, &new_width, &new_height, &new_channels, up_thres, low_thres);
            break;
        }
        case 9:{
            double cx = width / 2.0;
            double cy = height / 2.0;
            
            if(affine_transform_type == 1){
                 // 1. scaling ( 1.5x on  X and 0.8x on Y)
                result = affine_transform(img, width, height, channels, 
                    1.5, 0, 0, 0.8, 0, 0, &new_width, &new_height, affine_transform_type);
            }else if(affine_transform_type == 2){
                // 2. Rotation
                double angle = 180.0 * M_PI / 180.0;
                double cos_theta = cos(angle);
                double sin_theta = sin(angle);
                double tx = cx - cos_theta * cx + sin_theta * cy;
                double ty = cy - sin_theta * cx - cos_theta * cy;

                result = affine_transform(img, width, height, channels,
                    cos_theta, -sin_theta, sin_theta, cos_theta, tx, ty, &new_width, &new_height, affine_transform_type);
            }else if(affine_transform_type == 3){
                // 3. Translation (100 pixels to right 100 pixels down)
                result = affine_transform(img, width, height, channels, 
                    1, 0, 0, 1, 100, 50, &new_width, &new_height, affine_transform_type);
                
            }else if(affine_transform_type == 4){
                // 4. Shearing (X, Y)
                result = affine_transform(img, width, height, channels,
                    1, 0.5, 0.3, 1, 0, 0, &new_width, &new_height, affine_transform_type);
            }
        
            

            break;
        }
        default:
            printf("Incorrect operation choice.\n");
            stbi_image_free(img);
            return 1;
    }

    // 7. Results saving
    printf("Enter output file name (with extension: .jpg, .png, .bmp): ");    
    if (scanf("%255s", output_path) != 1) {
        fprintf(stderr, "Invalid filename input.\n");
        return 1;
    }


    char *ext = get_file_extension(output_path);
    for (int i = 0; ext[i]; i++) ext[i] = tolower(ext[i]);  

    int saved = 0;

    if (strcmp(ext, "jpg") == 0 || strcmp(ext, "jpeg") == 0) {
        saved = stbi_write_jpg(output_path, new_width, new_height, new_channels, result, 100);
    } else if (strcmp(ext, "png") == 0) {
        saved = stbi_write_png(output_path, new_width, new_height, new_channels, result, new_width * new_channels);
    } else if (strcmp(ext, "bmp") == 0) {
        saved = stbi_write_bmp(output_path, new_width, new_height, new_channels, result);
    } else {
        printf("Unsupported file extension. Supported: jpg, png, bmp\n");
        stbi_image_free(img);
        free(result);
        return 1;
    }

    if (saved) {
        printf("Image saved as %s\n", output_path);
    } else {
        fprintf(stderr, "Failed to save image.\n");
    }
    // 8. Free
    stbi_image_free(img);
    free(result);

    return 0;
}
