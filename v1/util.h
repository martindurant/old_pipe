#include "analyze_io.h"
/** @file
 * \brief Main collection of general-purpose utility routines.
 **/
typedef struct {
		/** \brief Number of columns in image */
		int cols;
		/** \brief Number of rows in image */
		int rows;
		/** \brief Number of slices in image */
		int slices;
		/** \brief Array of voxel intensities
		 *
		 * Voxel (x,y,z) is at position  z*rows*cols + y*cols + x.
		 */
		float *data;
		/** \brief Distance between centres of adjacent voxels */
		float sepx, sepy, sepz;
		/** \brief Minimum and maximum intensity values (from file) */
		float minval,maxval;
} image;

typedef struct {
		 /** \brief Nx6 array storing input (x,y,z) and output (x,y,z)
		  * coordinates for each node.
		  */
        float **nodemap;
		 /** \brief N-dimensional array storing intensity shift for each node.
		  */
        float *A;
		 /** \brief Uncertain, but required by image_register program. */
        float *MAP;
		 /** \brief Pixel separation between adjacent nodes. */
        int gridsize;
		 /**\brief dimension of node grid. */
        int xdim, ydim, zdim;
		 /** \brief dimension of source image. */
        int gridx, gridy, gridz;
		 /** \brief Total number of nodes. */
        int numnodes;
} map;


image *alloc_image(int cols,int rows,int slices);
void free_image(image *delimage);

image *interpolate_slices(image *in_image,float out_sep);
