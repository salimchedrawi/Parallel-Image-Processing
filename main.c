/*

 Parallel Image Processing

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void datread (char *filename, void *vx, int nx, int ny);
void pgmwrite(char *filename, void *vx, int nx, int ny);

#define M 192
#define N 360

#define Px 2
#define Py 2

#define P (Px*Py)

#define MP M/Px
#define NP N/Py

#define MAXITER   1500
#define PRINTFREQ  200

#define CHECKFREQ 100
#define THRESHOLD 0.1

int main (int argc, char **argv)
{
  float old[MP+2][NP+2], new[MP+2][NP+2], edge[MP+2][NP+2];

  float masterbuf[M][N];
  float buf[MP][NP];
  float buf2[M][N];

  int i, j, iter, maxiter,globali,globalj,location;
  char *filename;

  int rank, size, next, prev;
  int count=0;
  
  float PixelValue=0.0;
  float gdelta=0.0;
  float delta=0.0;
  float gdeltaall=0.0;
  float avgPixelValue=0.0;
  double start=0.0;
  double end=0.0;
  double timer=0.0;



  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  
  MPI_Comm cart_comm;
  
  MPI_Status status;
  MPI_Request req;

  if(size != P)
    {
     if (rank == 0) printf("ERROR: size = %d, P = %d\n", size, P);
      MPI_Finalize();
      exit(-1);
    }

  
/**********Decomposition***********/
    
    start = MPI_Wtime();
      
      MPI_Datatype blockType;
      MPI_Type_vector(M/Px,N/Py,N,MPI_FLOAT,&blockType);
      MPI_Type_commit(&blockType);



  int dims[2],periods[2],coords[2];

  dims[0]=Px; 
  dims[1]= Py;
  periods[0]=periods[1]=0;

  int ierr;
  int ierr2;

  ierr = MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1, &cart_comm);

  if (ierr != MPI_SUCCESS) printf("Cart create error!\n");


  MPI_Cart_coords(cart_comm, rank, 2, coords);



  int down_rank,up_rank,right_rank,left_rank;
  MPI_Cart_shift(cart_comm,0,1,&right_rank,&left_rank);
  MPI_Cart_shift(cart_comm,1,1,&down_rank,&up_rank);

  
  end = MPI_Wtime();
  
  timer = end - start;
  
  printf("Decompose: %f\n",timer);
  
  start = MPI_Wtime();
  
  /**********Read***********/

  if(rank == 0)
    {
      printf("Processing %d x %d image on %d processes\n", M, N, P);
      printf("Number of iterations = %d\n", MAXITER);

      filename = "edge192x360.dat";
      printf("\nReading <%s>\n", filename);
      datread(filename, masterbuf, M, N);
      printf("\n");
    }

  end = MPI_Wtime();
  
  timer = end - start;
  
  printf("Read: %f\n",timer);
  
  /**********Distribute***********/
  
    start = MPI_Wtime();

	ierr2 = MPI_Bcast(masterbuf, M*N, MPI_FLOAT, 0, MPI_COMM_WORLD);
       if (ierr2 != MPI_SUCCESS) printf("Send broadcast error!\n");



  for(i=0;i<MP;i++){
    for(j=0;j<NP;j++){

      globali = i + coords[0] * M/Px;
      globalj = j + coords[1] * N/Py;


      buf[i][j] = masterbuf[globali][globalj];

	  }
	  }

  end = MPI_Wtime();
  
  timer = end - start;
  
  printf("Distribute: %f\n",timer);
  
 /**********Array Swaps***********/
  
  start = MPI_Wtime();

  for (i=0;i<MP+2;i++)
    {
      for (j=0;j<NP+2;j++)
	{
	  edge[i][j]=0.0;
	}
    }
      
  for (i=1;i<MP+1;i++)
    {
      for (j=1;j<NP+1;j++)
	{
	  edge[i][j]=buf[i-1][j-1];
	}
    }

  for (i=0;i<MP+2;i++)
    {
      for (j=0;j<NP+2;j++)
	{
	  old[i][j]=edge[i][j];
	}
    }
    
    end = MPI_Wtime();
    
     timer = end - start;
  
  printf("Array Swaps %f\n",timer);


  MPI_Datatype haloType;
  MPI_Type_vector(MP,1,NP+2,MPI_FLOAT,&haloType);
  MPI_Type_commit(&haloType);
  
  
 /**********Reconstruction***********/
  
  start = MPI_Wtime();

  for (iter=1;iter<=MAXITER; iter++)
    {
      if(iter%PRINTFREQ==0)
	{
	  if(rank==0)
	    {
	      printf("Iteration %d\n", iter);
	    }
	}

  /**********Non-Blocking Communication for Boundary Elements***********/

		
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
		
		




   for (i=1;i<MP+1;i++)
     {
       for (j=1;j<NP+1;j++)
         {
           new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]
			      - edge[i][j]);
	 }
     }
     
  /**********Calculation of average pixel Value***********/     
     
     for(i=1;i<=MP;i++){
		for(j=1;j<=NP;j++)
		{
			delta += ((new[i][j] - old[i][j]) * (new[i][j] - old[i][j]));
			
			PixelValue += new[i][j];
				
		}
		}
		
	PixelValue /=(MP*NP);
	   
	MPI_Reduce(&PixelValue,&avgPixelValue,1,MPI_FLOAT,MPI_SUM,0,cart_comm);

	if(rank==0){
		if(iter%PRINTFREQ == 0)
		{
		avgPixelValue /=size;
		printf("Iteration number: %d, and average value of pixels: %f\n",iter,avgPixelValue);
		}
		}
			
						 		
	   
     
	
   for (i=1;i<MP+1;i++)
     {
       for (j=1;j<NP+1;j++)
         {
           old[i][j]=new[i][j];
	 }
     }
     
   /**********Image Accuracy Check***********/    
     
      MPI_Allreduce(&delta, &gdeltaall, 1, MPI_FLOAT, MPI_SUM, cart_comm);
	    
	   
	    gdelta = sqrt(gdeltaall/(M*N));
	    
     
     count++;
     
      if(gdelta < THRESHOLD)
	    {
	     printf("Iteration: %d\n\n", iter);
	     break;
	    }
     
     
     	     
	    PixelValue=0.0;
	    
	    gdeltaall=0.0;
	    
	    delta=0.0;
	    
	    gdelta=0.0;
	    
	    avgPixelValue=0.0;
	    
	    }
	    end = MPI_Wtime();
	    
	     timer = end - start;
  
  	printf("Reconstruction: %f\n",timer);
	
	timer = timer/count;
	
	printf("Average time per iteration: %f\n", timer);

   for (i=1;i<MP+1;i++)
     {
       for (j=1;j<NP+1;j++)
         {
           buf[i-1][j-1]=old[i][j];
	 }
     }
    

/**********Write***********/

 start = MPI_Wtime();

for(i=0;i<MP;i++){
    for(j=0;j<NP;j++){

      globali = i + coords[0] * M/Px;
      globalj = j + coords[1] * N/Py;


      buf2[globali][globalj] = buf[i][j];

	  }
	  }



 MPI_Allreduce(buf2,masterbuf,M*N,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);


  if (rank==0)
    {
      printf("\nFinished %d iterations\n", iter-1);
    }


  if (rank == 0)
    {
      filename="image192x360.pgm";
      printf("\nWriting <%s>\n", filename); 
      pgmwrite(filename, masterbuf, M, N);
    }
    
    end = MPI_Wtime();
    
     timer = end - start;
  
  printf("Write: %f\n",timer);

  MPI_Finalize();
} 
