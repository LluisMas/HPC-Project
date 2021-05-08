#!/bin/bash

## Specifies the interpreting shell for the job.
#$ -S /bin/bash

## Specifies that all environment variables active within the qsub utility be exported to the context of the job.
#$ -V

## Specifies the parallel environment
#$ -pe smp 2

## Execute the job from the current working directory.
#$ -cwd

## The  name  of  the  job.
#$ -N ConvolutionOMP_static

##Passes an environment variable to the job
#$ -v  OMP_NUM_THREADS=2

## In this line you have to write the command that will execute your application.

gcc convolution_omp_row_test.c -o 2_static_omp_convolution -fopenmp
for j in "25x25_random" "3x3_Edge" "49x49_random" "5x5_Sharpen" "99x99_random"
do
echo "Kernel $j.txt"
  echo "$IMAGES_PATH/im04.ppm $KERNEL_PATH/kernel${j}.txt > omp_img04_kernel${j}.txt"
  ./2_static_omp_convolution $IMAGES_PATH/im04.ppm $KERNEL_PATH/kernel${j}.txt /state/partition1/out.ppm > 2_static_omp_img04_kernel${j}.txt
#./2_static_omp_convolution $IMAGES_PATH/im04.ppm $KERNEL_PATH/kernel25x25_random.txt /state/partition1/out.ppm > test.txt
done
