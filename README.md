# CVlibC

**CVlibC** is a lightweight image processing library written in C that implements basic computer vision algorithms such as convolution, Gaussian blur, median filtering, grayscale conversion, Sobel operator, and Canny edge detection.

## 📦 Contents

- `CVlibC.h` — the main library with image processing functions realization
- `user_friendly_main.c` — the main library file containing image processing functions
- `Makefile` *(опционально)* —фва для компиляции проекта (можно добавить)

## 🚀 Functions

- Image convolution with custom kernels
- Gaussian blur
- Median filter
- Converting an image to grayscale
- Canny edge detection 
- Various Affine transformations
- Various kernels for convolution

## 🔧 Dependencies

Image I/O is handled via `stb_image` and `stb_image_write` — single-header libraries that are already included:

- [`stb_image.h`](https://github.com/nothings/stb/blob/master/stb_image.h)
    
- [`stb_image_write.h`](https://github.com/nothings/stb/blob/master/stb_image_write.h)
    

Download these files and place them in the same directory as the source code if not already present.

### Compilation

To build the project manually, run:

```
gcc user_friendly_main.c -o main -lm
```

