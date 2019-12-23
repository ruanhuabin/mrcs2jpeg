/*******************************************************************
 *       Filename:  mrcs2jpeg.c                                     
 *                                                                 
 *    Description:                                        
 *                                                                 
 *        Version:  1.0                                            
 *        Created:  11/07/2019 01:56:31 PM                                 
 *       Revision:  none                                           
 *       Compiler:  gcc                                           
 *                                                                 
 *         Author:  Ruan Huabin                                      
 *          Email:  ruanhuabin@tsinghua.edu.cn                                        
 *        Company:  Dep. of CS, Tsinghua Unversity                                      
 *                                                                 
 *******************************************************************/

/*
 * <setjmp.h> is used for the optional error recovery mechanism shown in
 * the second part of the example.
 */
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "jpeglib.h"
typedef struct _mrc_header__
{
    /* Number of columns, rows, and sections */
    int nx;
    int ny;
    int nz;

    /* Type of value in image:
     * 0 =unsigned or signed bytes depending on flag in imodStamp,only unsigned bytes before IMOD 4.2.23
     * 1 = signed short integers (16 bits)
     * 2 = float
     * 3 = short * 2, (used for complex data)
     * 4 = float * 2, (used for complex data)
     * 6 = unsigned 16-bit integers (non-standard)
     * 16 = unsigned char * 3 (for rgb data, non-standard)
     * */
    int mod;

    /* Start point of sub image */
    int nxstart;
    int nystart;
    int nzstart;

    /* Grid size in X, Y, and Z */
    int mx;
    int my;
    int mz;

    /* Cell size: pixel spacing = xlen/max, ylen/my, zlen/mz */
    float xlen;
    float ylen;
    float zlen;

    /* cell angles */
    float alpha;
    float beta;
    float gamma;

    int mapc;
    int mapr;
    int maps;

    float amin;
    float amax;
    float amean;

    int ispg;
    int next;

    char otherparts[1024 - 24 * 4];

} mrc_header_t;

typedef struct _mrc_fmt_
{
    mrc_header_t header;
    float *imageData;
} mrc_fmt_t;

void print_mrc_header(mrc_header_t *header)
{
    printf("Number of colums  : nx  = %d\n", header->nx);
    printf("Number of rows    : ny  = %d\n", header->ny);
    printf("Number of sections: nz  = %d\n", header->nz);
    printf("image data type   : mod = %d\n", header->mod);
}

void getMeanStd(float *data, size_t n, float &mean, float &std)
{
    float sum = 0.0f;
    for(size_t i = 0; i < n; i ++)
    {
        sum += data[i];
    }

    mean = sum / n;

    float variance = 0.0f;
    for(size_t i = 0; i < n; i ++)
    {
        variance += pow((data[i] - mean), 2);
    }

    variance /= n;

    std = sqrt(variance);
}


void preProcessImage(float *data, size_t n, float pMean, float pStd)
{
    float mean;
    float std;
    getMeanStd(data, n, mean, std);

    float maxValue = mean + (pMean + pStd) * std;
    float minValue = mean + (pMean -pStd) * std;
    for(size_t i = 0; i < n; i ++)
    {
        if(data[i] > maxValue)
        {
            printf("max value is changed: orig = %f, new = %f, max = %f\n", data[i], maxValue, maxValue);
            data[i] = maxValue;
        }

        if(data[i] < minValue)
        {
            printf("min value is changed: orig = %f, new = %f, min = %f\n", data[i], minValue, minValue);
            data[i] = minValue;
        }
    }
}
void getMinMax(float *data, int n, float& min, float& max, int& minIndex, int& maxIndex)
{
    min = data[0];
    max = data[0];
    minIndex = 0;
    maxIndex = 0;
    for(int i = 1; i < n; i ++)
    {
        if(data[i] < min)
        {
            min = data[i];
            minIndex = i;
        }

        if(data[i] > max)
        {
            max = data[i];
            maxIndex = i;
        }
    }
}


void addByAbsMin(float *data, int n, float& min)
{
    min = fabsf(min);
    for(int i = 0; i < n; i ++)
    {
        data[i] += min;
    }
}

void divByMax(float *data, int n, float&max)
{
    for(int i = 0; i < n; i ++)
    {
        data[i] /= max;
    }

}

void scaleToUncharByMul255(float *data, int n, unsigned char *out)
{
    for(int i = 0; i < n; i ++)
    {
        out[i] = (unsigned char)(data[i] * 255);
    }
}

void printFloatData(float *data, int start, int end, const char *hint)
{
    printf("%s:  ", hint);
    for(int i = start; i < end; i ++)
    {
        printf("%f ", data[i]);
    }

    printf("\n");
}

void printUncharData(unsigned char *data, int start, int end, const char *hint)
{
    printf("%s:  ", hint);
    for(int i = start; i < end; i ++)
    {
        printf("%d ", data[i]);
    }

    printf("\n");
}

/*
 * Sample routine for JPEG compression.  We assume that the target file name
 * and a compression quality factor are passed in.
 */
void write_JPEG_file (char * filename, int quality, unsigned char *image_buffer, int image_height, int image_width)
{
  /* This struct contains the JPEG compression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   * It is possible to have several such structures, representing multiple
   * compression/decompression processes, in existence at once.  We refer
   * to any one struct (and its associated working data) as a "JPEG object".
   */
  struct jpeg_compress_struct cinfo;
  /* This struct represents a JPEG error handler.  It is declared separately
   * because applications often want to supply a specialized error handler
   * (see the second half of this file for an example).  But here we just
   * take the easy way out and use the standard error handler, which will
   * print a message on stderr and call exit() if compression fails.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct jpeg_error_mgr jerr;
  /* More stuff */
  FILE * outfile;		/* target file */
  JSAMPROW row_pointer[1];	/* pointer to JSAMPLE row[s] */
  int row_stride;		/* physical row width in image buffer */

  /* Step 1: allocate and initialize JPEG compression object */

  /* We have to set up the error handler first, in case the initialization
   * step fails.  (Unlikely, but it could happen if you are out of memory.)
   * This routine fills in the contents of struct jerr, and returns jerr's
   * address which we place into the link field in cinfo.
   */
  cinfo.err = jpeg_std_error(&jerr);
  /* Now we can initialize the JPEG compression object. */
  jpeg_create_compress(&cinfo);

  /* Step 2: specify data destination (eg, a file) */
  /* Note: steps 2 and 3 can be done in either order. */

  /* Here we use the library-supplied code to send compressed data to a
   * stdio stream.  You can also write your own code to do something else.
   * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
   * requires it in order to write binary files.
   */
  if ((outfile = fopen(filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    exit(1);
  }
  jpeg_stdio_dest(&cinfo, outfile);

  /* Step 3: set parameters for compression */

  /* First we supply a description of the input image.
   * Four fields of the cinfo struct must be filled in:
   */
  cinfo.image_width = image_width; 	/* image width and height, in pixels */
  cinfo.image_height = image_height;
  cinfo.input_components = 1;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_GRAYSCALE; 	/* colorspace of input image */
  /* Now use the library's routine to set default compression parameters.
   * (You must set at least cinfo.in_color_space before calling this,
   * since the defaults depend on the source color space.)
   */
  jpeg_set_defaults(&cinfo);
  /* Now you can set any non-default parameters you wish to.
   * Here we just illustrate the use of quality (quantization table) scaling:
   */
  jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

  /* Step 4: Start compressor */

  /* TRUE ensures that we will write a complete interchange-JPEG file.
   * Pass TRUE unless you are very sure of what you're doing.
   */
  jpeg_start_compress(&cinfo, TRUE);

  /* Step 5: while (scan lines remain to be written) */
  /*           jpeg_write_scanlines(...); */

  /* Here we use the library's state variable cinfo.next_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   * To keep things simple, we pass one scanline per call; you can pass
   * more if you wish, though.
   */
  row_stride = image_width * 1;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    row_pointer[0] = & image_buffer[cinfo.next_scanline * row_stride];
    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }

  /* Step 6: Finish compression */

  jpeg_finish_compress(&cinfo);
  /* After finish_compress, we can close the output file. */
  fclose(outfile);

  /* Step 7: release JPEG compression object */

  /* This is an important step since it will release a good deal of memory. */
  jpeg_destroy_compress(&cinfo);

  /* And we're done! */
}

int main ( int argc, char *argv[] )
{ 

    if (argc != 6)
    {
        printf("Usage: %s <mrc file name> <output directory> <pMean> <pStd> <sNormalized>\n", argv[0]);
        printf("e.g: %s image.mrcs output 0 3 1\n", argv[0]);
        exit(-1);
    }

    FILE *mrcFile = fopen(argv[1], "rb");

    float pMean = atof(argv[3]);
    float pStd = atof(argv[4]);
    int sNormalized = atoi(argv[5]);

    if (mrcFile == NULL)
    {
        printf("Open input file [%s] failed: [%s:%d]", argv[1], __FILE__, __LINE__);
        exit(-1);
    }


    mrc_fmt_t mrcInstance;
    fread(&mrcInstance.header, sizeof(mrc_header_t), 1, mrcFile);
    print_mrc_header(&mrcInstance.header);
    int nx = mrcInstance.header.nx;
    int ny = mrcInstance.header.ny;
    int nz = mrcInstance.header.nz;
    int mod = mrcInstance.header.mod;

    printf("[nx, ny, nz] = [%d, %d, %d]\n", nx, ny, nz);

    if(mod != 2)
    {
        fprintf(stderr, "Error: Data type in this mrc file is not float, please check.\n");
        exit(1);
    }

    size_t n = nx * ny;
    mrcInstance.imageData = (float *)malloc(n * sizeof(float));
    unsigned char *grayData = (unsigned char *)malloc(n * sizeof(unsigned char));
    float min = 0.0f;
    float max = 0.0f;
    int minIndex = 0;
    int maxIndex = 0;
    char cmd[1024];
    strcpy(cmd, "mkdir -p ");
    strcat(cmd, argv[2]);
    system(cmd);
    
    char outputJpegFileName[256];

    for(int i = 0; i < nz; i ++)
    {
        fread(mrcInstance.imageData, sizeof(float), nx * ny, mrcFile);

        if(sNormalized != 0)
        {
            preProcessImage(mrcInstance.imageData, n, pMean, pStd);
        }

        getMinMax(mrcInstance.imageData, n, min, max, minIndex, maxIndex); 
        printf("[min, max, minIndex, maxIndex] = [%f, %f, %d, %d]\n", min, max, minIndex, maxIndex);
        if(min < 0)
        {
            fprintf(stderr, "Warn: min value is < 0, all the elements will be added by value of min\n");
            /*printFloatData(mrcInstance.imageData, 0, 10, "before add by min");*/
            addByAbsMin(mrcInstance.imageData, n, min);
            /*printFloatData(mrcInstance.imageData, 0, 10, "after add by min");*/
            max = max + fabsf(min);
            min = 0.0f;
            fprintf(stderr, "[newMin, newMax]= [%f, %f]\n", min, max);

            if(max < 1e-6)
            {
                fprintf(stderr, "Error: The max value of input file is zero, please check.\n");
                exit(1);
            }

        }

        /*printFloatData(mrcInstance.imageData, 0, 10, "before div by max");*/
        divByMax(mrcInstance.imageData, n, max);
        /*printFloatData(mrcInstance.imageData, 0, 10, "after div by max");*/
        scaleToUncharByMul255(mrcInstance.imageData, n, grayData);
        /*printUncharData(grayData, 0, 10, "after scale to 255");*/

        sprintf(outputJpegFileName, "%s/%d.jpeg", argv[2], (i+1));
        write_JPEG_file(outputJpegFileName, 100, grayData, nx, ny);

    }
    

    free(mrcInstance.imageData);
    free(grayData);
    fclose(mrcFile);

    return EXIT_SUCCESS;
}

