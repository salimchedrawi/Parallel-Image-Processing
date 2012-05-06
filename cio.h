/*

cio.h
is written for the Message-Passing Programming Coursework
MSc in HPC, EPCC, 2010
Author: Salim CHEDRAWI

*/
#ifndef CIO_H
#define CIO_H


/*
 *  Routine to read an "edges" data file into a 2D floating point array
 *  x[nx][ny]. Because of the way C handles (or fails to handle!)
 *  multi-dimensional arrays we have to cast the pointer to void.
 */

void datread(char *filename, void *vx, int nx, int ny);



/*
 *  Routine to write a PGM image file from a 2D floating point array
 *  x[nx][ny]. Because of the way C handles (or fails to handle!)
 *  multi-dimensional arrays we have to cast the pointer to void.
 */

void pgmwrite(char *filename, void *vx, int nx, int ny);

#endif /* CIO_H */
