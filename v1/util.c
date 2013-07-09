#include "util.h"
#include <stdlib.h>
#include <string.h>
#include "analyze_io.h"
#include "splines.h"
#include <math.h>
#include <stdio.h>


/**
 * \brief Allocate memory for image structure
 *
 * Free with free_image
 **/
image *alloc_image(int cols, int rows, int slices) {
		image *retimage;

		retimage = (image *)malloc(sizeof(image));
        retimage->cols = cols;
	    retimage->rows = rows;
	    retimage->slices = slices;
	    retimage->data = (float*)malloc(sizeof(float)*cols*rows*slices);

		return retimage;
}

/**
 * \brief Safely deallocate memory for image structure
 **/
void free_image(image *delimage) {
		free(delimage->data);
		free(delimage);
}

/**
 * \brief Applies a cubic spline interpolation to the slices of the input image.
 *
 * Produces a slice separation given by out_sep.
 **/
image *interpolate_slices(image *in_image, float out_sep) {
		image *out_image;
		double in_sep;
		int i, j, k;
		float *fx;
		float *fy;
		float *fy2;
		float y[1];


		if (fabs(in_image->sepx - in_image->sepy) > 0.01) {
			printf("Slices in input file don't have square pixels.\n");
			return NULL;
		}
		in_sep = in_image->sepz;

		out_image = alloc_image(in_image->cols,in_image->rows,(int)(in_image->slices-1)*(in_sep/out_sep)+1);
		fx = malloc(sizeof(float)*in_image->slices);
		fy = malloc(sizeof(float)*in_image->slices);
		fy2 = malloc(sizeof(float)*in_image->slices);

		for (k = 0; k<in_image->slices; k++)
				fx[k] = k*in_sep;

		for (i = 0; i<in_image->rows; i++) {
			for (j = 0; j<in_image->cols; j++) {
					for (k = 0; k < in_image->slices; k++)
							fy[k] = in_image->data[k*in_image->rows*in_image->cols + i*in_image->cols + j];
					spline(fx, fy, in_image->slices, 0, 0, fy2);
					for (k = 0; k < out_image->slices; k++) {
							splint(fx,fy,fy2,in_image->slices,k*out_sep,y);
							out_image->data[k*out_image->rows*out_image->cols + i*out_image->cols + j] = y[0];
					}
			}
		}
		out_image->sepx = in_image->sepx;
		out_image->sepy = in_image->sepy;
		out_image->sepz = out_sep;
		free(fx);
		free(fy);
		free(fy2);

		return out_image;
}

