#include "correlation.h"
#include <array>
#include <iostream>

/************************ TODO 2 **************************/
/*
*	INPUT:
*		origImg:		the original image,
*		imgWidth:		the width of the image
*		imgHeight:		the height of the image
*						the image is arranged such that
*						origImg[3*(row*imgWidth+column)+0],
*						origImg[3*(row*imgWidth+column)+1],
*						origImg[3*(row*imgWidth+column)+2]
*						are R, G, B values for pixel at (column, row).
*
*      kernel:			the 2D filter kernel,
*		knlWidth:		the width of the kernel
*		knlHeight:		the height of the kernel
*
*		scale, offset:  after correlating the kernel with the origImg,
*						each pixel should be divided by scale and then added by offset
*
*		selection:      a byte array of the same size as the image,
*						indicating where in the original image should be filtered, e.g.,
*						selection[k] == 1 ==> pixel k should be filtered
*                      selection[k] == 0 ==> pixel k should NOT be filtered
*                      a special case is selection is a NULL pointer, which means all the pixels should be filtered.
*
*  OUTPUT:
*		rsltImg:		the filtered image of the same size as original image.
*						it is a valid pointer ( allocated already ).
*/

void image_filter(double* rsltImg, const unsigned char* origImg, const unsigned char* selection,
				  int imgWidth, int imgHeight,
				  const double* kernel, int knlWidth, int knlHeight,
				  double scale, double offset)
{
	double rsltpxl[3];

	if (selection != NULL) //If non-NULL, must check every pixel
	{
		for (int y = 0; y < imgHeight; y++) {
			for (int x = 0; x < imgWidth; x++) {
				if (selection[y * imgWidth + x] == 1) {
					pixel_filter(rsltpxl, x, y, origImg, imgWidth, imgHeight, kernel, knlWidth, knlHeight, scale, offset);

					rsltImg[3 * (y * imgWidth + x) + 0] = rsltpxl[0];
					rsltImg[3 * (y * imgWidth + x) + 1] = rsltpxl[1];
					rsltImg[3 * (y * imgWidth + x) + 2] = rsltpxl[2];
				}
				else {
					rsltImg[3 * (y * imgWidth + x) + 0] = origImg[3 * (y * imgWidth + x) + 0];
					rsltImg[3 * (y * imgWidth + x) + 1] = origImg[3 * (y * imgWidth + x) + 1];
					rsltImg[3 * (y * imgWidth + x) + 2] = origImg[3 * (y * imgWidth + x) + 2];
				}
			}
		}
	}
	else //If NULL, can move through the whole image.
	{
		for (int y = 0; y < imgHeight; y++) {
			for (int x = 0; x < imgWidth; x++) {
				pixel_filter(rsltpxl, x, y, origImg, imgWidth, imgHeight, kernel, knlWidth, knlHeight, scale, offset);

				rsltImg[3 * (y * imgWidth + x) + 0] = rsltpxl[0];
				rsltImg[3 * (y * imgWidth + x) + 1] = rsltpxl[1];
				rsltImg[3 * (y * imgWidth + x) + 2] = rsltpxl[2];
			}
		}
	}
}

/************************ END OF TODO 2 **************************/


/************************ TODO 3 **************************/
/*
*	INPUT:
*      x:				a column index,
*      y:				a row index,
*		origImg:		the original image,
*		imgWidth:		the width of the image
*		imgHeight:		the height of the image
*						the image is arranged such that
*						origImg[3*(row*imgWidth+column)+0],
*						origImg[3*(row*imgWidth+column)+1],
*						origImg[3*(row*imgWidth+column)+2]
*						are R, G, B values for pixel at (column, row).
*
*      kernel:			the 2D filter kernel,
*		knlWidth:		the width of the kernel
*		knlHeight:		the height of the kernel
*
*		scale, offset:  after correlating the kernel with the origImg,
*						the result pixel should be divided by scale and then added by offset
*
*  OUTPUT:
*		rsltPixel[0], rsltPixel[1], rsltPixel[2]:
*						the filtered pixel R, G, B values at row y , column x;
*/

void pixel_filter(double rsltPixel[3], int x, int y, const unsigned char* origImg, int imgWidth, int imgHeight,
				  const double* kernel, int knlWidth, int knlHeight,
				  double scale, double offset)
{
	int colors = 3;
	for (int color = 0; color < 3; ++color) {
		double result = 0.0000;
		for (int kernelX = 0; kernelX < knlWidth; ++kernelX) {
			int windowX = x - knlWidth/2 + kernelX;
			for (int kernelY = 0; kernelY < knlHeight; ++kernelY) {
				int windowY = y - knlHeight/2 + kernelY;
				if ((0 <= windowX && windowX < imgWidth) && (0 <= windowY && windowY < imgHeight)) {
					result += (kernel[kernelY*knlWidth+kernelX]) * (origImg[3*(windowY*imgWidth+windowX)+color]);
				}
			}
		}
		rsltPixel[color] = result/scale + offset;
	}
}

/************************ END OF TODO 3 **************************/