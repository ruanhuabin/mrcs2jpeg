/*******************************************************************
 *       Filename:  cDisp_3DMap.c                                     
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
#include <getopt.h>
#include "jpeglib.h"
#include "readsection.h"
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
            data[i] = maxValue;
        }

        if(data[i] < minValue)
        {
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
void divByMax2(float *data, int n, float &min, float &max)
{
    float diff = max - min;
    for(int i = 0; i < n; i ++)
    {
        data[i] = (data[i] - min) / diff;
    }
}

void scaleToUncharByMul255(float *data, int n, unsigned char *out)
{
    int i = 0;
    for(i = 0; i < n; i ++)
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

void usage(const char *program)
{
    printf("Usage:\n\t%s --input <input file name> --output <output directory> --pMean <float value> --pStd <float value> --pAxis <x|y|z> --pSlice <int num> [sNormalized] [sMRC|sJPEG|SRAW]\n", program);
    printf("\n\te.g: %s --input image.mrcs --output output --pMean 0.0 --pStd 3.0 --pAxis z --pSlice 0 --sNormalized --sJPEG\n", program);
}
int main( int argc, char *argv[] )
{ 

    if(argc < 4)
    {
        usage(argv[0]);
        return 1;
    }
    int c;
    float pMean = 0.0f;
    float pStd = 1.0f;
    int sNormalized = 0;
    int sJPEG = 0;
    int sRAW = 0;
    int sMRC = 0;
    int pSlice = 0;
    char pAxis = 'z';
    char inputFilename[256];
    char outputDir[256];
    char outputFilePrefix[256];
    memset(inputFilename, '\0', sizeof(inputFilename));
    memset(outputDir, '\0', sizeof(outputDir));
    memset(outputFilePrefix, '\0', sizeof(outputFilePrefix));
    while(1)
    {
        static struct option long_options[] =
        {
          {"sNormalized" , no_argument       , 0 , 'n'} ,
          {"sMRC"        , no_argument       , 0 , 'm'} ,
          {"sJPEG"       , no_argument       , 0 , 'j'} ,
          {"sRAW"        , no_argument       , 0 , 'r'} ,
          {"pMean"       , required_argument , 0 , 'a'} ,
          {"pStd"        , required_argument , 0 , 's'} ,
          {"pAxis"       , required_argument , 0 , 'x'} ,
          {"pSlice"      , required_argument , 0 , 'y'} ,
          {"input"       , required_argument , 0 , 'i'} ,
          {"output"      , required_argument , 0 , 'o'} ,
          {"prefix"      , required_argument , 0 , 'p'} ,
          {"help"        , no_argument       , 0 , 'h'} ,
          {0             , 0                 , 0 , 0 }
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "mnjrhp:a:s:i:o:x:y:", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c)
        {
        case 'n':
          sNormalized = 1;
          break;

        case 'm':
          sMRC = 1;
          break;

        case 'j':
          sJPEG = 1;
          break;
        case 'r':
          sRAW = 1;
          break;

        case 'a':
          pMean = atof(optarg);
          break;

        case 's':
          pStd = atof(optarg);
          break;
        case 'x':
          /*pAxis= atoi(optarg);*/
          pAxis = optarg[0];
          break;
        case 'y':
          pSlice = atoi(optarg);
          break;
        case 'i':
          /*printf("input file name: %s\n", optarg);*/
          strcpy(inputFilename, optarg); 
          break;

        case 'o':
          strcpy(outputDir, optarg); 
          break;

        case 'p':
          strcpy(outputFilePrefix, optarg); 
          break;

        case 'h':
          /* getopt_long already printed an error message. */
          usage(argv[0]);
          return 1;

        default:
          usage(argv[0]);
          return 1;
        }

    
    }

    /*printf("input = %s, output = %s, sNormalized = %d, sJPEG = %d, sRAW = %d, sMRC = %d,  pMean = %f, pStd = %f, pAxis = %d, pSlice = %d\n", inputFilename, outputDir, sNormalized, sJPEG, sRAW, sMRC, pMean, pStd, pAxis, pSlice);*/

    /*return 0;*/

    
    char cmd[1024];
    strcpy(cmd, "mkdir -p ");
    strcat(cmd, outputDir);
    system(cmd);

    /*printf("create output dir success: %s\n", cmd);*/

    float *im;
    int dimx;
    int dimy;
    MRCHeader mrcHeader;
    MRC mrc;
    int ret = mrc.open(inputFilename, "rb");
    if(ret == 0)
    {
        fprintf(stderr, "Open file [%s] failed\n", inputFilename);
        return 1;
    }
    mrc.getHeader(&mrcHeader);
    readMRCSection2(&im, dimx, dimy, mrc, pAxis, pSlice);

    /**
     *  write nx, ny, nz to meta.txt
     */

    char metaFilename[1024];
    /*strcat(metaFilename, "/meta.txt");*/
    sprintf(metaFilename, "%s/%s_meta.txt", outputDir, outputFilePrefix);
    /*printf("meta file name: %s\n", metaFilename);*/
    FILE *metaFile = fopen(metaFilename, "w");
    char pDimXStr[64];
    char pDimYStr[64];
    char pDimZStr[64];

    if(sJPEG == 1)
    {
        fputs("pDataType=sUChar\n", metaFile);
    }
    else
    {
        fputs("pDataType=sReal32\n", metaFile);
    }

    sprintf(pDimXStr, "pDimX=%d\n", dimx);
    sprintf(pDimYStr, "pDimY=%d\n", dimy);
    sprintf(pDimZStr, "pDimZ=%d\n", 1);
    
    fputs(pDimXStr, metaFile);
    fputs(pDimYStr, metaFile);
    fputs(pDimZStr, metaFile);

    fclose(metaFile);

    size_t n = dimx * dimy;
    unsigned char *grayData = (unsigned char *)malloc(n * sizeof(unsigned char));
    float min = 0.0f;
    float max = 0.0f;
    int minIndex = 0;
    int maxIndex = 0;
    
    char outputImageFileName[256];
    float mean = 0.0f;
    float std = 0.0f;

    if(sNormalized == 1)
    {
        //for(size_t i = 0; i < nz; i ++)
        //{
        //fread(mrcInstance.imageData, sizeof(float), n, mrcFile);
        //getMeanStd(mrcInstance.imageData, n, mean, std);
        getMeanStd(im, n, mean, std);
        for(size_t j = 0; j < n; j ++)
        {
            //mrcInstance.imageData[j] = (mrcInstance.imageData[j] - mean) / std * pStd + pMean;
            im[j] = (im[j] - mean) / std * pStd + pMean;
        }

        //Write current image as mrc file format
        if(sMRC == 1)
        {
        
            sprintf(outputImageFileName, "%s/%s_%d.mrc", outputDir, outputFilePrefix, pSlice);
            FILE *f = fopen(outputImageFileName, "w");

            mrcHeader.nx = dimx;
            mrcHeader.ny = dimy;
            mrcHeader.nz = 1;
            fwrite(&mrcHeader, sizeof(MRCHeader),1, f);
            fwrite(im, sizeof(float), n, f);
            fclose(f); 
        }
        else if(sRAW == 1)
        {
            sprintf(outputImageFileName, "%s/%s_%d.raw", outputDir, outputFilePrefix, pSlice);
            FILE *f = fopen(outputImageFileName, "w");
            fwrite(im, sizeof(float), n, f);
            fclose(f); 
        }
        //}

    }
    else if(sJPEG == 1)
    {
        float maxValue = 0.0f;
        float minValue = 0.0f;
        getMeanStd(im, n, mean, std);
        maxValue = mean + std * (pMean + pStd);
        minValue = mean + std * (pMean - pStd);
        for(size_t j = 0; j < n; j ++)
        {
            if(im[j] > maxValue)
            {
                im[j] = maxValue;
            }

            if(im[j] < minValue)
            {
                im[j] = minValue;
            }

        }
        getMinMax(im, n, min, max, minIndex, maxIndex); 
        if(min < 0)
        {
            /*fprintf(stderr, "Warn: min value is < 0, all the elements will be added by value of min\n");*/
            addByAbsMin(im, n, min);
            max = max + fabsf(min);
            min = 0.0f;
            if(max < 1e-6)
            {
                fprintf(stderr, "Warn: The max value of input file is zero, please check.\n");
                /*exit(1);*/
            }
        }

        divByMax2(im, n, min, max);
        scaleToUncharByMul255(im, n, grayData);
        
        sprintf(outputImageFileName, "%s/%s_%d.jpeg", outputDir, outputFilePrefix, pSlice);
        write_JPEG_file(outputImageFileName, 100, grayData, dimy, dimx);
    }
    free(grayData);
    mrc.close();
    delete[] im;
    
    return EXIT_SUCCESS;
}

