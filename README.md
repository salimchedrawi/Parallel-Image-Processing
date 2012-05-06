Parallel-Image-Processing
=========================

The main purpose of this project is to use image processing to reconstruct an image using a data file that was the output of an edge detection algorithm such as Jacobi or Gauss-Seidel. However the aim is to do two-dimensional lattice-based calculation using a two-dimensional domain decomposition that uses non-blocking communications.


1 Introduction
=========================


The program has been written on C and works correctly. Many tests have been done in order to experiment, study, and evaluate the performance of the program. 
As a result, the speed (time of execution) and performance (avg. time per iteration) of the software has increased with the increase of the number of processors. 


2 Experimental design and Implementation
=========================

2.1 Method
=========================

The program inputs a data file with information that contains edge-detection data that is outputted from an algorithm. The edge-detection algorithm takes a greyscale image of a supercomputer which is of size M x N and gives the data of the edges. This data is then used by the program to reconstruct the image and outputs the image in a PGM format. 
First there should be a decomposition for the number of processors P. A two dimensional decomposition is chosen. Px is the number of processors in the first dimension. Py is the number of processors in the second dimension. As such, P = Px x Py.
Since the image size is of size M x N, then the data of the edges should be read and put in memory and then M/Px x N/Py is dispersed to each processor.
The program then runs calculations and uses non-blocking communication for halo swapping and reconstructs the image through various iterations
The result is then outputted as a PGM file.

2.2 Decomposition
=========================

There should be a two dimensional decomposition for P in order to decrease communication due to large overheads. Therefore P is split into Px and Py. Px is the number of processors in the first dimension. Py is the number of processors in the second dimension. As such, P = Px x Py. 
Since the size of the image is M x N, then the challenge is to find Px and Py so as to minimize the communications overhead found in halo swapping. Therefore, for practicality, in order to minimize (Px * N) + (Py * M), M should be divisibly by Px and N should be divisble by Py.

2.3 Reading and Distributing Data
=========================

During the first approach for the 1-D decomposition, the data file was read and distributed using MPI_Scatter since each processor should work on a segment of the data in the edges file. However this cannot be done for the 2-D decomposition. I found two good methods to implement this. One method is where each process reads the edges data file and takes the segment that it needs. The other method, which I chose in my program, is where the master process (rank 0) reads the edges data file and puts it into an array. It then broadcasts the entire array to each processor. 
In order to broadcast the data, the master process sends it using MPI_Bcast. Since it is a broadcast then we do not need a MPI_Recv, so no synchronization is needed.
In order for the master process to read the data file, the function 'datread' was used. However in order to give the required segment of data to each processor respectively, new respective coordinates where used which segments the data by M/Px and M/Py according to the respective cartesian plane of each processor. The pseudo for such a step is the following:

loop over i=0,MP and j=0,NP
globali = i + coords[0] * M/Px;
globalj = j + coords[1] * N/Py;
buf[i][j] = masterbuf[globali][globalj];

2.4 Image Reconstruction
=========================

In order to reconstruct the image, four static arrays are created: float buf of size [MP][NP]. This is used to store the required data of the edges of the respective processor. The arrays float old, float new, and float edge that are of size [MP+2][NP+2]. These are used to reconstruct the image. We notice that they are bigger than the array buf, this is because halo values for swaps are stored in them. As well as the array float masterbuf which is of size [M][N]. This array is used to put in the resulted data after calculations and create the PGM image file for display.
Since the program uses a 2-D decomposition, in order to do calculations to finally reconstruct the image, halo swapping should be used. As such, for each process to know its neighbours, a 2-D Cartesian topology is used with of course non-periodic boundary conditions. 
Since halo swaps are involved at the edges, then in order to send and receive the rows and columns the cartesian shifting must be used. For instance, to send and receive the rows in the 'old' array, the up_rank and down_rank must be used to shift up and down between ranks. As for the columns, left_rank and right_rank should be used to shift right and left between ranks. 
A derived datatype is used in order to send and receive the columns for the halo swaps between ranks as such: 

MPI_Datatype haloType;
MPI_Type_vector(MP,1,NP+2,MPI_FLOAT,&haloType);
MPI_Type_commit(&haloType);

2.4.1 Iterative Process of Image Reconstruction
=========================
The iterative process goes through a number of MAXITER iterations. The more number of iterations, the better the image is reconstructed.
More computation than communication is used during the iterative process so as to re-construct the image.
The iterative process goes as follows:

2.4.1.1	Boundary Elements
=========================
1) This is a non-blocking communication procedure where the boundary elements are swapped. The boundary element swaps are finalized and ensured using MPI_Wait. The code is the following:

MPI_Issend(&old[1][NP], 1,haloType,up_rank, 0,cart_comm, &req);
MPI_Recv(&old[1][0], 1,haloType, down_rank, 0,cart_comm, &status);
MPI_Wait(&req, &status);
MPI_Issend(&old[1][1], 1,haloType,down_rank, 0,cart_comm, &req);
MPI_Recv(&old[1][NP+1], 1,haloType, up_rank, 0,cart_comm, &status);
MPI_Wait(&req, &status);
MPI_Issend(&old[MP][1], NP,MPI_FLOAT,left_rank, 0,cart_comm, &req);
MPI_Recv(&old[0][1], NP,MPI_FLOAT, right_rank, 0,cart_comm, &status);
MPI_Wait(&req, &status);
MPI_Issend(&old[1][1], NP,MPI_FLOAT,right_rank, 0,cart_comm, &req);
MPI_Recv(&old[MP+1][1], NP,MPI_FLOAT, left_rank, 0,cart_comm, &status);
MPI_Wait(&req, &status);

2) The non boundary and boundary elements are calculated. The code is the following:

loop over i=1,MP and j=1,NP
new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]- edge[i][j]);

2.4.1.2	Copying of arrays
=========================

3) In order to the store the calculation results, the old array is copied to the new array without the halos. The code is the following:

loop over i=1,MP+1 and j=1,NP+1
old[i][j]=new[i][j];

2.4.1.3	Calculation of the average pixel values and image accuracy check
=========================

4) The average pixel values of the image that was reconstructed are calculated. They are then printed at a specific frequency. The code is the following:

loop over i=1,MP and j=1,NP
PixelValue += new[i][j];
PixelValue /=(MP*NP);
MPI_Reduce(&PixelValue,&avgPixelValue,1,MPI_FLOAT,MPI_SUM,0,cart_comm);
if(rank==0){
if(iter%PRINTFREQ == 0)
{
avgPixelValue /=size;
printf("Iteration number: %d, and average value of pixels: %f\n",iter,avgPixelValue);
}
		}

5) We check the accuracy of the image by comparing it to a specific THRESHOLD. If it is not accurate enough, the loop breaks.
Here, the ? parameter is calculated. It is then compared to a THRESHOLD to see the accuracy of the image. If it is not that accurate, the loop is terminated. This is put in the loop to minimize overheads.

loop over i=1,MP and j=1,NP
delta+= ((new[i][j] - old[i][j]) * (new[i][j] - old[i][j]));
MPI_Allreduce(&delta, &gdeltaall, 1, MPI_FLOAT, MPI_SUM, cart_comm);
gdelta = sqrt(gdeltaall/(M*N);
if(gdelta < THRESHOLD)
{
printf("Iteration: %d\n\n", iter);
break;  }

2.4.2 Copy of data after loop termination
=========================

When the loop is completed, the data (without the halos) in the old array is copied to the buf array. 

loop over i=1,MP+1 and j=1,NP+1
buf[i-1][j-1]=old[i][j];

2.5 Final image data
=========================

In order to construct the final image, all the data generated by each processor should be aggregated at the masterbuf array. As such, since every data segment of each processor are put in the buf array, each processor will put its buff array to a temp array buff2 of size [M][N] using the same new respective coordinates that were used to segment the image in the beginning by M/Px and M/Py according to the respective cartesian plane of each processor. As such, the array will contain within its specific space, the required data for each processor. The rest of the spaces will be zero.
However in order to give the required segment of data to each processor respectively, new respective coordinates where used which segments the data by M/Px and M/Py according to the respective cartesian plane of each processor. The pseudo for such a step is the following:

loop over i=0,MP and j=0,NP
globali = i + coords[0] * M/Px;
globalj = j + coords[1] * N/Py;
 buf2[globali][globalj] = buf[i][j];

Then in order to join all the temporary arrays together and have one big masterbuf array with all the information of all the processors at the master process an MPI_Allreduce is used.
At the end, at the master process, the masterbuf is passed to pgmwrite in order to write out the final reconstructed image.
