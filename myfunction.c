#include <stdbool.h>
#include <stdlib.h>
// AVOIDING CALCULATIONS FOR KNOWN SIZES AND FREQUENT CALCULATIONS BY SETTING AS CONSTS.
#define KERNEL_SIZE                     3
#define HALF_KERNEL_SIZE                1 //KERNEL_SIZE/2
//AVOID UNNECESSARY CALLS FOR A FUNC WITH OVER 40,000,000 CALLS LOL
#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))


typedef struct {
   unsigned char red;
   unsigned char green;
   unsigned char blue;
} pixel;

typedef struct {
    int red;
    int green;
    int blue;
    // int num;
} pixel_sum;


int calcIndex(int i, int j, int n) {
	return ((i)*(n)+(j));
}

/*
 * initialize_pixel_sum - Initializes all fields of sum to 0
 * NO NEEDED AFTER OPTIMIZATION.
 */
/*void initialize_pixel_sum(pixel_sum *sum) {
	sum->red = sum->green = sum->blue = 0;
	// sum->num = 0;
	return;
}*/

/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) {

	// divide by kernel's weight
	sum.red /= kernelScale;
	sum.green = sum.green / kernelScale;
	sum.blue = sum.blue / kernelScale;

	// truncate each pixel's color values to match the range [0,255]
	current_pixel->red = (unsigned char) (min(max(sum.red, 0), 255));
	current_pixel->green = (unsigned char) (min(max(sum.green, 0), 255));
	current_pixel->blue = (unsigned char) (min(max(sum.blue, 0), 255));
	return;
}

/*
* sum_pixels_by_weight - Sums pixel values, scaled by given weight
*/
static void sum_pixels_by_weight(pixel_sum *sum, pixel p, int weight) {

	sum->red += ((int) p.red) * weight;
	sum->green += ((int) p.green) * weight;
	sum->blue += ((int) p.blue) * weight;

	// sum->num++;
	return;
}

/*
 *  Applies kernel for pixel at (i,j)
 */
static pixel applyKernel (int dim, int i, int j, pixel *src, int kernel[KERNEL_SIZE][KERNEL_SIZE], int kernelScale, bool filter) {

    int kRow, kCol, ii, jj;
    //INITIALIZING STRUCT USING INIT LIST. MINOR EFFECT
	pixel_sum sum = {0, 0, 0};
	pixel current_pixel;
	int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
	int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
	int min_row, min_col, max_row, max_col;
	pixel loop_pixel;

	//initialize_pixel_sum(&sum);
    //INSTEAD OF DOING THE SAME CALC ALL OVER AND OVER
    int iLoopCounter = min(i+1, dim-1);
    int jLoopCounter = min(j+1, dim-1);
    //AVOIDED SCANNING MATRIX TWICE BY COPING LOOP'S CONTEXT INTO THE NEXT LOOP
    //AND AVOIDING UNNECESSARY SCANS.
    if(!filter) {
        for(ii = max(i-1, 0); ii <= iLoopCounter; ++ii) {
            for(jj = max(j-1, 0); jj <= jLoopCounter; ++jj) {

                // compute row index in kernel
                if (ii < i) {
                    kRow = 0;
                } else if (ii > i) {
                    kRow = 2;
                } else {
                    kRow = 1;
                }

                // compute column index in kernel
                if (jj < j) {
                    kCol = 0;
                } else if (jj > j) {
                    kCol = 2;
                } else {
                    kCol = 1;
                }

                // apply kernel on pixel at [ii,jj]
                sum_pixels_by_weight(&sum, src[calcIndex(ii, jj, dim)], kernel[kRow][kCol]);
            }
        }
	}

/*    iLoopCounter = min(i+1, dim-1);
    jLoopCounter = min(j+1, dim-1);*/

	else {
		// find min and max coordinates
		for(ii = max(i-1, 0); ii <= iLoopCounter; ii++) {
			for(jj = max(j-1, 0); jj <= jLoopCounter; jj++) {
                // compute row index in kernel
                if (ii < i) {
                    kRow = 0;
                } else if (ii > i) {
                    kRow = 2;
                } else {
                    kRow = 1;
                }

                // compute column index in kernel
                if (jj < j) {
                    kCol = 0;
                } else if (jj > j) {
                    kCol = 2;
                } else {
                    kCol = 1;
                }

                // apply kernel on pixel at [ii,jj]
                sum_pixels_by_weight(&sum, src[calcIndex(ii, jj, dim)], kernel[kRow][kCol]);
				// check if smaller than min or higher than max and update
				loop_pixel = src[calcIndex(ii, jj, dim)];
				if ((((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue)) <= min_intensity) {
					min_intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
					min_row = ii;
					min_col = jj;
				}
				if ((((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue)) > max_intensity) {
					max_intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
					max_row = ii;
					max_col = jj;
				}
			}
		}
		// filter out min and max
		sum_pixels_by_weight(&sum, src[calcIndex(min_row, min_col, dim)], -1);
		sum_pixels_by_weight(&sum, src[calcIndex(max_row, max_col, dim)], -1);
	}

	// assign kernel's result to pixel at [i,j]
	assign_sum_to_pixel(&current_pixel, sum, kernelScale);
	return current_pixel;
}

/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
void smooth(int dim, pixel *src, pixel *dst, int kernel[KERNEL_SIZE][KERNEL_SIZE], int kernelScale, bool filter) {

	int i, j;
	for (i = HALF_KERNEL_SIZE ; i < dim - HALF_KERNEL_SIZE; i++) {
		for (j =  HALF_KERNEL_SIZE ; j < dim - HALF_KERNEL_SIZE ; j++) {
			dst[calcIndex(i, j, dim)] = applyKernel(dim, i, j, src, kernel, kernelScale, filter);
		}
	}
}

void charsToPixels(Image *charsImg, pixel* pixels) {

	int row, col;
	/*for (row = 0 ; row < m ; row++) {
		for (col = 0 ; col < n ; col++) {

			pixels[row*n + col].red = image->data[3*row*n + 3*col];
			pixels[row*n + col].green = image->data[3*row*n + 3*col + 1];
			pixels[row*n + col].blue = image->data[3*row*n + 3*col + 2];
		}
	}*/
    // CONVERTING:
    // pixels[row*n + col].red = image->data[3*row*n + 3*col];
    // = pixels[row*n + col].red = image->data[3*(row*n + col)];
    // TOTAL: RUNS N * M TIMES.
    // BASICALLY IT GOES THROUGH EACH INT IN RANGE 0 TO M * N
    // LET i BE (row * n + col)
    // SUCH THAT:
    // pixels[row*n + col].red = image->data[3*row*n + 3*col]
    // = pixels[i].red = image->data[3*(i)];
    int iterations = m * n, curPos = 0;
    char *dataPtr = image -> data;
    for (int i = 0; i < iterations; i++) {
        pixels[i].red = *dataPtr;
        ++dataPtr;
        pixels[i].green = *dataPtr;
        ++dataPtr;
        pixels[i].blue = *dataPtr;
        ++dataPtr;
    }

}

void pixelsToChars(pixel* pixels, Image *charsImg) {

	int row, col;
	/*for (row = 0 ; row < m ; ++row) {
		for (col = 0 ; col < n ; ++col) {
			image->data[3*row*n + 3*col] = pixels[row*n + col].red;
			image->data[3*row*n + 3*col + 1] = pixels[row*n + col].green;
			image->data[3*row*n + 3*col + 2] = pixels[row*n + col].blue;
		}
	}*/
    // SIMILAR OPTIMIZATION TO CHARS_TO_PIXELS(...)
    int iterations = m * n, curPos = 0;
    char *dataPtr = image -> data;
    for (int i = 0; i < iterations; i++) {
        *dataPtr = pixels[i].red;
        ++dataPtr;
        *dataPtr = pixels[i].green;
        ++dataPtr;
        *dataPtr = pixels[i].blue;
        ++dataPtr;
    }
}

void copyPixels(pixel* src, pixel* dst) {

	int row, col, add;
	for (row = 0 ; row < m ; row++) {
		for (col = 0 ; col < n ; col++) {
            // CALC ROW*N+COL ONCE, INSTEAD OF 6 TIMES. SEEMS AS MINOR UPDATE. (1)
            add = row * n + col;
			dst[add].red = src[add].red;
			dst[add].green = src[add].green;
			dst[add].blue = src[add].blue;
		}
	}
}

void doConvolution(Image *image, int kernel[KERNEL_SIZE][KERNEL_SIZE], int kernelScale, bool filter) {

	pixel* pixelsImg = malloc(m*n*sizeof(pixel));
	pixel* backupOrg = malloc(m*n*sizeof(pixel));

	charsToPixels(image, pixelsImg);
	copyPixels(pixelsImg, backupOrg);

	smooth(m, backupOrg, pixelsImg, kernel, kernelScale, filter);

	pixelsToChars(pixelsImg, image);

	free(pixelsImg);
	free(backupOrg);
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

	/*
	* [1, 1, 1]
	* [1, 1, 1]
	* [1, 1, 1]
	*/
	int blurKernel[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

	/*
	* [-1, -1, -1]
	* [-1, 9, -1]
	* [-1, -1, -1]
	*/
	int sharpKernel[3][3] = {{-1,-1,-1},{-1,9,-1},{-1,-1,-1}};

	if (flag == '1') {	
		// blur image
		doConvolution(image, blurKernel, 9, false);

		// write result image to file
		writeBMP(image, srcImgpName, blurRsltImgName);	

		// sharpen the resulting image
		doConvolution(image, sharpKernel, 1, false);
		
		// write result image to file
		writeBMP(image, srcImgpName, sharpRsltImgName);	
	} else {
		// apply extermum filtered kernel to blur image
		doConvolution(image, blurKernel, 7, true);

		// write result image to file
		writeBMP(image, srcImgpName, filteredBlurRsltImgName);

		// sharpen the resulting image
		doConvolution(image, sharpKernel, 1, false);

		// write result image to file
		writeBMP(image, srcImgpName, filteredSharpRsltImgName);	
	}
}

