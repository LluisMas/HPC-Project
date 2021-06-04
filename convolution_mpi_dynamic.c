//
//  convolution.c
//
//
//  Created by Josep Lluis Lerida on 11/03/15.
//  Updated by Vitor da Silva on 17/04/2021
//
// This program allows you to apply the convolution to an image file with a * .ppm extension.
// The program receives the file with the source image, the file with the kernel for the convolution and the path of the output file.
// The 2D matrix of the image is represented by a 1D vector for each R, G and B channel. The convolution is applied for each channel separately.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define CHUNK 10

// Structure for storing the content of an image.
struct imagenppm{
    int height;
    int width;
    char *comment;
    int maxcolor;
    int P;
    int *R;
    int *G;
    int *B;
};
typedef struct imagenppm* ImagenData;

// Structure for storing the contents of a kernel.
struct structkernel{
    int kernelX;
    int kernelY;
    float *vkern;
};
typedef struct structkernel* kernelData;

// Definition of the main functions: reading an image, duplicating the image, reading the kernel, convolution, saving an image in PPM format, etc.
ImagenData readImage(char* name);
kernelData readKernel(char* name);
ImagenData duplicateImageData(ImagenData src);
int convolve2D(int*, int*, int, int, float*, int, int);
int saveFile(ImagenData Img, char* name);

// This Function allows us to read a ppm file and place the information: encoding, size, RGB, etc. of the image in a structure.
// The value of each RGB pixel in the 2D image is saved and represented by 1D vectors for each channel.
ImagenData  readImage(char* name){
    FILE *fp;
    char c;
    char comment[300];
    int i=0;
    ImagenData Img=NULL;

    // The ppm file is opened
    fp=fopen(name,"r");
    if(!fp){
        perror("Error");
    }
    else{
        // We reserve memory for the Image structure
        Img=(ImagenData) malloc(sizeof(struct imagenppm));

        //The magic number is read and saved
        fscanf(fp,"%c%d ",&c,&(Img->P));
        // The comment is read and saved character by character
        while((c=fgetc(fp))!= '\n'){comment[i]=c;i++;}
        Img->comment = calloc(strlen(comment),sizeof(char));
        strcpy(Img->comment,comment);
        // Read and save the width, height and maximum color
        fscanf(fp,"%d %d %d",&Img->width,&Img->height,&Img->maxcolor);
        // Memory is reserved in R, G and B according to width and height
        // And the values of R, G and B of the file are assigned
        if ((Img->R=calloc(Img->width*Img->height,sizeof(int))) == NULL) {return NULL;}
        if ((Img->G=calloc(Img->width*Img->height,sizeof(int))) == NULL) {return NULL;}
        if ((Img->B=calloc(Img->width*Img->height,sizeof(int))) == NULL) {return NULL;}
        for(i=0;i<Img->width*Img->height;i++){
            fscanf(fp,"%d %d %d ",&Img->R[i],&Img->G[i],&Img->B[i]);
        }
        fclose(fp);
    }
    return Img;
}

// This function allows us to read a kernel from a file, the kernel is represented by a 1D vector.
kernelData readKernel(char* name){
    FILE *fp;
    int i=0;
    kernelData kern=NULL;

    //Opening ppm file
    fp=fopen(name,"r");
    if(!fp){
        perror("Error: ");
    }
    else{
        //We reserve memory for the structure that the kernel will store
        kern=(kernelData) malloc(sizeof(struct structkernel));

        // We read the dimensions of the kernel.
        fscanf(fp,"%d,%d,", &kern->kernelX, &kern->kernelY);
        kern->vkern = (float *)malloc(kern->kernelX*kern->kernelY*sizeof(float));

        //kernel reading
        for (i=0;i<(kern->kernelX*kern->kernelY)-1;i++){
            fscanf(fp,"%f,",&kern->vkern[i]);
        }
        fscanf(fp,"%f",&kern->vkern[i]);
        fclose(fp);
    }
    return kern;
}


// This function allows you to copy the main data from the original image data structure into a second structure.
ImagenData duplicateImageData(ImagenData src){
    unsigned int imageX, imageY;

    // We reserve memory for the target Image structure
    ImagenData dst=(ImagenData) malloc(sizeof(struct imagenppm));
    //Magic number is copied
    dst->P=src->P;
    // Comment is copied
    dst->comment = calloc(strlen(src->comment),sizeof(char));
    strcpy(dst->comment,src->comment);
    // Width, height and maximum color are copied
    imageX=src->width;
    imageY=src->height;
    dst->width=imageX;
    dst->height=imageY;
    dst->maxcolor=src->maxcolor;
    // Memory is reserved in R, G and B according to width and height
    if ((dst->R=calloc(imageX*imageY,sizeof(int))) == NULL) {return NULL;}
    if ((dst->G=calloc(imageX*imageY,sizeof(int))) == NULL) {return NULL;}
    if ((dst->B=calloc(imageX*imageY,sizeof(int))) == NULL) {return NULL;}
    memcpy(dst->R, src->R, (imageX*imageY)*sizeof(int));
    memcpy(dst->G, src->G, (imageX*imageY)*sizeof(int));
    memcpy(dst->B, src->B, (imageX*imageY)*sizeof(int));
    return dst;
}

// This function stores the new Image data in a ppm file.
int saveFile(ImagenData  Img, char* name){
    int i;
    FILE *fp;
    // Resulting image is created
    if (!(fp=fopen(name,"w"))) {
        printf("Error opening the result file: %s\n",name);
        return -1;
    }
    //The magic number, comment, width, height and maximum color are written to the file
    fprintf(fp,"P%d\n%s\n%d %d\n%d\n",Img->P,Img->comment,Img->width,Img->height,Img->maxcolor);
    //Pixels are written
    for(i=0;i<Img->width*Img->height;i++){
            fprintf(fp,"%d %d %d ",Img->R[i],Img->G[i],Img->B[i]);
            if (i%Img->height==0) fprintf(fp,"\n");
    }
    fclose(fp);
    return 0;
}


///////////////////////////////////////////////////////////////////////////////
// 2D convolution
// 2D data are usually stored in computer memory as contiguous 1D array.
// So, we are using 1D array for 2D data.
// 2D convolution assumes the kernel is center originated, which means, if
// kernel size 3 then, k[-1], k[0], k[1]. The middle of index is always 0.
// The following programming logics are somewhat complicated because of using
// pointer indexing in order to minimize the number of multiplications.
//
//
// signed integer (32bit) version:
///////////////////////////////////////////////////////////////////////////////
int convolve2D(int* in, int* out, int dataSizeX, int dataSizeY,
                float* kernel, int kernelSizeX, int kernelSizeY)
{
    int i, j, m, n;
    int *inPtr, *inPtr2, *outPtr;
    float *kPtr;
    int kCenterX, kCenterY;
    int rowMin, rowMax;                             // to check boundary of input array
    int colMin, colMax;                             //
    float sum;                                      // temp accumulation buffer

    // check validity of params
    if(!in || !out || !kernel) return -1;
    if(dataSizeX <= 0 || kernelSizeX <= 0) return -1;

    // find center position of kernel (half of kernel size)
    kCenterX = (int)kernelSizeX / 2;
    kCenterY = (int)kernelSizeY / 2;

    // init working  pointers
    inPtr = inPtr2 = &in[dataSizeX * kCenterY + kCenterX];  // note that  it is shifted (kCenterX, kCenterY),
    outPtr = out;
    kPtr = kernel;

    for(i= 0; i < dataSizeY; i ++)                   // number of rows
    {
        // compute the range of convolution, the current row of kernel should be between these
        rowMax = i + kCenterY;
        rowMin = i - dataSizeY + kCenterY;

        for(j = 0; j < dataSizeX; ++j)              // number of columns
        {
            //printf("rowMax: %d, rowMin %d \n", rowMax, rowMin);
            // compute the range of convolution, the current column of kernel should be between these
            colMax = j + kCenterX;
            colMin = j - dataSizeX + kCenterX;
            sum = 0;                                // set to 0 before accumulate

            // flip the kernel and traverse all the kernel values
            // multiply each kernel value with underlying input data
            for(m = 0; m < kernelSizeY; ++m)        // kernel rows
            {
                // check if the index is out of bound of input array
                if(m <= rowMax && m > rowMin)
                {
                    for(n = 0; n < kernelSizeX; ++n)
                    {
                        // check the boundary of array
                        if(n <= colMax && n > colMin)
                            sum += *(inPtr - n) * *kPtr;

                        ++kPtr;                     // next kernel
                    }
                }
                else
                    kPtr += kernelSizeX;            // out of bound, move to next row of kernel

                inPtr -= dataSizeX;                 // move input data 1 raw up
            }

            // convert integer number
            if(sum >= 0) *outPtr = (int)(sum + 0.5f);
            //else *outPtr = (int)(sum - 0.5f)*(-1);
            // For using with image editors like GIMP or others...
            else *outPtr = (int)(sum - 0.5f);
            // For using with a text editor that read ppm images like libreoffice or others...
            // else *outPtr = 0;

            kPtr = kernel;                          // reset kernel to (0,0)
            inPtr = ++inPtr2;                       // next input
            ++outPtr;                               // next output
        }
      }

    return 0;
}

int main(int argc, char **argv)
{
    int rank, size;
    int dimensions[2];
    char hostname[256];
    int namelen;
    ImagenData  source=NULL, output=NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(hostname, &namelen);   // get CPU name

    if(argc != 4)
    {
        printf("Usage: %s <image-file> <kernel-file> <result-file>\n", argv[0]);

        printf("\n\nError, missing parameters:\n");
        printf("format: image_file kernel_file result_file\n");
        printf("- image_file : source image path (*.ppm)\n");
        printf("- kernel_file: kernel path (text file with 1D kernel matrix)\n");
        printf("- result_file: result image path (*.ppm)\n\n");
        MPI_Finalize();
        return -1;
    }

    struct timeval tim;
    gettimeofday(&tim, NULL);
    double t1=tim.tv_sec+(tim.tv_usec/1000000.0);

    if (rank == 0) // Only the master process needs to read all of that
    {
      //Read the source image
      if ( (source=readImage(argv[1]))==NULL) {
          MPI_Finalize();
          return -1;
      }

      gettimeofday(&tim, NULL);
      double t2=tim.tv_sec+(tim.tv_usec/1000000.0);

      // Duplicate the image in a new structure that will contain the output image
      if ( (output=duplicateImageData(source)) == NULL) {
          MPI_Finalize();
          return -1;
      }

      dimensions[0] = source->height;
      dimensions[1] = source->width;

      printf("Height: %d, Width: %d \n", source->height, source->width);
    }

    MPI_Bcast(&dimensions, 2, MPI_INT, 0, MPI_COMM_WORLD);
    gettimeofday(&tim, NULL);
    double t3=tim.tv_sec+(tim.tv_usec/1000000.0);

    printf("[%d] Dimensions: H: %d W: %d, Size: %d\n", rank, dimensions[0], dimensions[1], size);
    printf("[%d] Going to read Kernel %s \n", rank, argv[2]);

    //Kernel reading
    kernelData kern=NULL;
    if ( (kern = readKernel(argv[2]))==NULL) {
        printf("[%d] Failed\n", rank);
        MPI_Finalize();
        return -1;
    }

    printf("[%d] Kernel done\n", rank);


    gettimeofday(&tim, NULL);
    double t4=tim.tv_sec+(tim.tv_usec/1000000.0);

    double starttime, endtime;
    if (rank == 0)
      starttime = MPI_Wtime();

    int remainingWork = dimensions[0];
    MPI_Status status;

    // Send initial work
    if (rank == 0)
    {
      int i;
      for (i = 0; i < size-1; i++)
      {
        int start = dimensions[0] - remainingWork;
        int amount = remainingWork >= CHUNK ? CHUNK : remainingWork;

        // Every process receives 3 vectors (R,B,G) of Width x CHUNK size
        int *sendArray = (int*)malloc(sizeof(int) * (dimensions[1]*CHUNK*3));
        MPI_Send(sendArray, dimensions[1]*CHUNK*3, MPI_INT, i+1, start+1, MPI_COMM_WORLD);

        remainingWork -= amount;
      }

      printf("Initial data sent \n");
    }

    int finishSignalsSent = 0;
    while (1)
    {
      if (rank == 0)
      {
        int flag = 0;
        int start = dimensions[0] - remainingWork;
        while (!flag)
        {
          if (remainingWork > 0)
          {
            // Do some work while waiting
            convolve2D(source->R + start, output->R + start, dimensions[1], 1, kern->vkern, kern->kernelX, kern->kernelY);
            convolve2D(source->G + start, output->G + start, dimensions[1], 1, kern->vkern, kern->kernelX, kern->kernelY);
            convolve2D(source->B + start, output->B + start, dimensions[1], 1, kern->vkern, kern->kernelX, kern->kernelY);
            start = dimensions[0] - --remainingWork;
          }
          // Check if workers are done
          MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status );
        }

        printf("Received TAG: %d\n", status.MPI_TAG);
        int *recv = (int*)malloc(sizeof(int) * (dimensions[1]*CHUNK*3));
        int currentSize = status.MPI_TAG > 0 ? dimensions[1]*CHUNK*3 : (status.MPI_TAG*-1);

        MPI_Recv(recv, currentSize, MPI_INT, MPI_ANY_SOURCE , status.MPI_TAG, MPI_COMM_WORLD, &status );

        int displacement  = status.MPI_TAG > 0 ? status.MPI_TAG : dimensions[1] - (status.MPI_TAG*-1);
        memcpy(output->R + displacement, recv, sizeof(int) * currentSize );
        int *sendArray = (int*)malloc(sizeof(int) * (dimensions[1]*CHUNK*3));
        if (remainingWork <= 0)
        {
          // Sending 0 to alert the worker that there is no more data
          MPI_Send(sendArray, dimensions[1]*CHUNK*3, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
          if (++finishSignalsSent == size-1)
          {
            printf("Done\n");
            break;
          }
        }
        else
        {
          int amount = remainingWork >= CHUNK ? CHUNK : remainingWork;
          int end = start + amount;
          remainingWork -= amount;

          MPI_Send(sendArray, dimensions[1]*CHUNK*3, MPI_INT, status.MPI_SOURCE, start+1, MPI_COMM_WORLD);
          printf("Start: %d End: %d Amount: %d Remaining %d \n", start, end, amount, remainingWork);
        }

      }
      else
      {
        int *recv = (int*)malloc(sizeof(int) * (dimensions[1]*CHUNK*3));
        int *sendBuff = (int*)malloc(sizeof(int) * (dimensions[1]*CHUNK*3));
        int i;

        MPI_Recv(recv, dimensions[1]*CHUNK*3, MPI_INT,MPI_ANY_SOURCE , MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        printf("[%d] Data reveiced \n", rank);
        if (status.MPI_TAG == 0)// Done
        {
          //printf("[%d] Done, bye\n", rank);
          break;
        }

        convolve2D(recv, sendBuff, dimensions[1], CHUNK, kern->vkern, kern->kernelX, kern->kernelY);
        convolve2D(recv+dimensions[1], sendBuff+dimensions[1], dimensions[1], CHUNK, kern->vkern, kern->kernelX, kern->kernelY);
        convolve2D(recv+dimensions[1]*2, sendBuff+dimensions[1] * 2, dimensions[1], CHUNK, kern->vkern, kern->kernelX, kern->kernelY);

        MPI_Send(sendBuff, dimensions[1]*CHUNK*3, MPI_INT, 0, status.MPI_TAG, MPI_COMM_WORLD);
      }
    }
    printf("[%d] Waiting on the barrier\n", rank );
    MPI_Barrier(MPI_COMM_WORLD);


/*
    convolve2D(outR, rec_bufR, dimensions[1], sendcount[rank], kern->vkern, kern->kernelX, kern->kernelY);
    convolve2D(outG, rec_bufG, dimensions[1], sendcount[rank], kern->vkern, kern->kernelX, kern->kernelY);
    convolve2D(outB, rec_bufB, dimensions[1], sendcount[rank], kern->vkern, kern->kernelX, kern->kernelY);

*/
    if (rank == 0)
    {
      endtime = MPI_Wtime();
      printf("Start: %f, \n", starttime);
      printf("End: %f, \n", endtime);
      printf("Elapsed: %f\n", endtime-starttime);

      // Image writing
      if (saveFile(output, argv[3])!=0) {
          printf("Error saving the image\n");
          MPI_Finalize();
          return -1;
      }
    }
    /*
    #pragma omp parallel num_threads(2)
    {
      #pragma omp sections
      {
        #pragma omp section
          convolve2D(source->R, output->R, source->width, source->height, kern->vkern, kern->kernelX, kern->kernelY);
        #pragma omp section
          convolve2D(source->G, output->G, source->width, source->height, kern->vkern, kern->kernelX, kern->kernelY);
        #pragma omp section
          convolve2D(source->B, output->B, source->width, source->height, kern->vkern, kern->kernelX, kern->kernelY);
      }
    }
*/
    gettimeofday(&tim, NULL);
    double t5=tim.tv_sec+(tim.tv_usec/1000000.0);


    gettimeofday(&tim, NULL);
    double t6=tim.tv_sec+(tim.tv_usec/1000000.0);
/*
    printf("Image: %s\n", argv[1]);
    printf("SizeX : %d\n", source->width);
    printf("SizeY : %d\n", source->height);
    printf("%.6lf seconds elapsed for Reading image file.\n", t2-t1);
    printf("%.6lf seconds elapsed for copying image structure.\n", t3-t2);
    printf("%.6lf seconds elapsed for Reading kernel matrix.\n", t4-t3);
    printf("%.6lf seconds elapsed for make the convolution.\n", t5-t4);
    printf("%.6lf seconds elapsed for writing the resulting image.\n", t6-t5);
    printf("%.6lf seconds elapsed\n", t6-t1);
*/
    MPI_Finalize();
    return 0;
}
