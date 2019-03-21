//Compile with 
//g++ -O3 -fopenmp mvmultiply_sparse.cpp -o mvmultiply_sparse

#include <iostream>
#include <cassert>
#include <omp.h>
#include <cmath>
#include <cstdlib>

//A struct to hold the data for a CRSMatrix
struct CRSMatrix
{
    //The values of the matrix elements, array length is nnz.
    double* value;
    //column[i] is the column index for the i:th element in values, length nnz.
    unsigned int* column;
    //RowPtr[i] is the first index that corresponds to row i. Length NRows+1. The last element is nnz.
    unsigned int* RowPtr;

    //Number of rows in the matrix
    unsigned int NRows;
    //Number of columns in the matrix
    unsigned int NCols;
    //Number of nonzero elements
    unsigned int nnz;

    //Allocate memory for the matrix
    CRSMatrix(unsigned int NRows,unsigned int NCols,unsigned int nnz):
        NRows(NRows),NCols(NCols),nnz(nnz)
    {
        value=new double[nnz];
        column=new unsigned int[nnz];
        RowPtr=new unsigned int[NRows+1];
    }

    //Free the allocated memory
    ~CRSMatrix()
    {
        delete[] value;
        value=NULL;
        delete[] column;
        column=NULL;
        delete[] RowPtr;
        RowPtr=NULL;
    }
};

//A struct to hold the data for a COOMatrix
struct COOMatrix
{
    //The values of the matrix elements, array length is nnz.
    double* value;
    //The column indices. column[i] is the column corresponding to value[i].
    unsigned int* column;
    //The row indices. row[i] is the row corresponding to value[i].
    unsigned int* row;

    //Number of rows in the matrix
    unsigned int NRows;
    //Number of columns in the matrix
    unsigned int NCols;
    //Number of nonzero elements
    unsigned int nnz;

    //Allocate memory for the matrix
    COOMatrix(unsigned int NRows,unsigned int NCols,unsigned int nnz):
	NRows(NRows),NCols(NCols),nnz(nnz)
    {
        value=new double[nnz];
        column=new unsigned int[nnz];
        row=new unsigned int[nnz];
    }

    //Free the allocated memory
    ~COOMatrix()
    {
        delete[] value;
        value=NULL;
        delete[] column;
        column=NULL;
        delete[] row;
        row=NULL;
    }
};


//Copy the data from a matrix given in the COO-format into another
//one given in the CRS-format.
void CopyCOOToCRS(COOMatrix & M1,CRSMatrix & M2)
{    
    assert(M1.NRows==M2.NRows && M1.NCols==M2.NCols && M1.nnz==M2.nnz);

    std::fill(M2.RowPtr,M2.RowPtr+M2.NRows+1,0);

    for(unsigned int i=0;i<M1.nnz;++i)
    {
        M2.RowPtr[M1.row[i]+1]++;
    }

    for(unsigned int i=1;i<=M1.NRows;++i)
    {
        M2.RowPtr[i]+=M2.RowPtr[i-1];
    }
    
    assert(M2.RowPtr[M2.NRows]==M1.nnz);

    for(unsigned int i=0;i<M1.nnz;++i)
    {
        unsigned int row=M1.row[i];
        M2.value[M2.RowPtr[row]]=M1.value[i];
        M2.column[M2.RowPtr[row]]=M1.column[i];
        M2.RowPtr[row]++;
    }
    
    for(unsigned int i=M2.NRows;i!=0;--i)
    {
        M2.RowPtr[i]=M2.RowPtr[i-1];
    }
    M2.RowPtr[0]=0;
    assert(M2.RowPtr[M2.NRows]==M2.nnz);
}

//Initialize a COO-matrix to contain the five-point stencil.
//We will use this matrix in our tests.
//You can access the data in the matrix structure using the construct matrix.field.
//For instance, if j is an integer, matrix.value[j], matrix.column[j], and matrix.row[j]
//are data for the jth element in the COO sparse matrix structure
void FillMatrix(COOMatrix & matrix)
{
  //j will be the counter for the matrix entries, increase it by one every time
  //you add a new entry
  unsigned int j = 0;
  unsigned int Nx = sqrt(matrix.NCols);
  unsigned int x = 0;
  unsigned int y = 0;
  //It is easiest to proceed row-by-row
  for (unsigned int row=0;row < matrix.NCols;++row)
    {
      //std::cout<<row<<std::endl;
      //x and y can be used to compute the grid position
      y = row/Nx;
      x = row - y*Nx;
      //Here comes the code for the entries of the matrix
        if(row < Nx){
            matrix.value[j] = 1;
            matrix.column[j] = row;
            matrix.row[j] = row;
            j++;
        }else if(row < matrix.NCols - Nx){
            if (row % Nx == 0 || row % Nx - 1 == 0){
                matrix.value[j] = 1;
                matrix.column[j] = x;
                matrix.row[j] = row;
                j++;
            }else{
                matrix.value[j] = -1;
                matrix.column[j] = x - Nx;
                matrix.row[j] = row;

                matrix.value[j + 1] = -1;
                matrix.column[j + 1] = x - 1;
                matrix.row[j + 1] = row;

                matrix.value[j + 2] = -4;
                matrix.column[j] = x;
                matrix.row[j] = row;

                matrix.value[j + 3] = -1;
                matrix.column[j] = x + 1;
                matrix.row[j] = row;

                matrix.value[j + 4] = -1;
                matrix.column[j] = x + Nx;
                matrix.row[j] = row;
                j += 5;
            }
        }else{
            matrix.value[j] = 1;
            matrix.column[j] = row;
            matrix.row[j] = row;
            j++;
        }
      }
  //Check that the number of entries if correct
  if (j != matrix.nnz) {
    std::cout<<"Error! COO matrix not properly filled!"<<std::endl;
    std::cout<<"Counting index and nnz differ."<<std::endl;
    std::cout<<"nnz: "<<matrix.nnz<<std::endl;
    std::cout<<"j: "<<j<<std::endl;
  }
}

//Print a matrix given in the COO format.
void PrintCOOMatrix(const COOMatrix & M)
{
    for(unsigned int i=0;i<M.nnz;++i)
    {
        std::cout<<M.row[i]<<" "<<M.column[i]<<" "<<M.value[i]<<std::endl;
    }
}

//Print a matrix given in the CRS format.
void PrintCRSMatrix(const CRSMatrix & M)
{
    for(unsigned int i=0;i<M.nnz;++i)
    {
        std::cout<<M.column[i]<<" ";
    }
    std::cout<<std::endl;

    for(unsigned int i=0;i<M.nnz;++i)
    {
        std::cout<<M.value[i]<<" ";
    }
    std::cout<<std::endl;

    for(unsigned int i=0;i<=M.NRows;++i)
    {
        std::cout<<M.RowPtr[i]<<" ";
    }
    std::cout<<std::endl;
}

//Put some test data into the array vector whose length is N.
void FillVector(unsigned int N,double* vector)
{
    for(unsigned int i=0;i<N;++i)
    {
      vector[i]=0.01*(rand() % 100);
    }
}

//Print the array vector whose length is N.
void PrintVector(unsigned int N,double* vector)
{
    for(unsigned int i=0;i<N;++i)
    {
        std::cout<<vector[i]<<" ";
    }
    std::cout<<std::endl;
}

//Return true if the vectors are the same within a tolerance, 
//and false otherwise.
bool CompareVectors(unsigned int N,double * v1,double * v2)
{
    const double tolerance=1e-10;
    for(unsigned int i=0;i<N;++i)
    {
        if(std::abs(v1[i]-v2[i])>tolerance)
        {
            return(false);
        }
    }

    return(true);
}

//Perform the multiplication matrix*vector and save the result in result.
//The length of the vectors must be the same as the number of rows in matrix.
void matvec_coo(const COOMatrix & matrix, double* vector, double* result)
{
    //Note that you can refer to the data contained in matrix using the
    //. operator, and index the arrays using the [] operator.
    //For example, the first element is matrix.value[0] and the corresponding
    //row and column are matrix.row[0] and matrix.column[0]. The number of rows
    //and columns in the matrix can be accessed as matrix.NRows and matrix.NCols
    //and the number of nonzero elements is matrix.nnz.
    std::fill(result,result+matrix.NRows,0);

    for (int k = 0; k < matrix.nnz; k = k + 1){
        result[matrix.row[k]] = result[matrix.row[k]] + matrix.value[k]*vector[matrix.column[k]];
    }
}

//Perform the multiplication matrix*vector and save the result in result.
//The length of the vectors must be the same as the number of rows in matrix.
void matvec_crs(const CRSMatrix & matrix,double* vector,double* result)
{   
    //omp_set_num_threads (2) ;
    for(unsigned int i=0;i<matrix.NRows;++i)
    {
        result[i]=0;
 
        for (int k = matrix.RowPtr[i]; k < matrix.RowPtr[i+1]; k = k + 1)
        {  
            result[i] = result[i] + matrix.value[k]*vector[matrix.column[k]];
      }  
   }  
    
    
}

void matvec_crs_par(const CRSMatrix & matrix,double* vector,double* result)
{   
    omp_set_num_threads (3) ;
    #pragma omp  parallel for
    for(unsigned int i=0;i<matrix.NRows;++i)
    {
        result[i]=0;
 
        for (int k = matrix.RowPtr[i]; k < matrix.RowPtr[i+1]; k = k + 1)
        {  
            result[i] = result[i] + matrix.value[k]*vector[matrix.column[k]];
      }  
   }  
    
    
}




int main()
{
    //Grid size in one dimension and corresponding spacing
    const unsigned int Nx=6000;
    const double h=1.0/(Nx-1);
    //Number of grid points = number of rows and columns in the matrix
    const unsigned int N=Nx*Nx;
    
    //Compute the number of nonzero elements for the five-points stencil
    //with Dirichlet boundary conditions and enter the result here instead of the
    //placeholder -1
    const unsigned int nnz= 5 * (Nx - 2) * (Nx - 2) + 2 * Nx + 2 * (Nx - 2); 
    
    //Create a matrix in the COO format and fill it with the stencil data.
    COOMatrix MCOO(N,N,nnz);
    FillMatrix(MCOO);

    //Create a CRS matrix and fill it with the same data as MCOO.
    CRSMatrix MCRS(N,N,nnz);
    CopyCOOToCRS(MCOO,MCRS);

    //Create a test vector to be multiplied by the matrices.
    double * vector=new double[N];
    FillVector(N,vector);

    //Allocate space for the results.
    double * result_coo=new double[N];
    double * result_crs=new double[N];
    
    
    
    int n_repeats = 5;
    //Do the multiplications and time.
    double time0=omp_get_wtime();
    
    //for(int i = 0; i < n_repeats;i++){ 
        matvec_crs(MCRS,vector,result_crs);
    //}

    double time1=omp_get_wtime();
    
    //for(int i = 0; i < n_repeats; i++){ 
        matvec_crs_par(MCRS,vector,result_crs);
    //}

    double time2=omp_get_wtime();

    matvec_coo(MCOO, vector, result_coo);
    //matvec_crs_par(MCRS,vector,result_crs);
    
    //Print results.
    std::cout<<"Time for Regular: "<<time1-time0<<std::endl;
    std::cout<<"Time for Parallel: "<<time2-time1<<std::endl;


    //Check that the rsults are the same.
    if(!CompareVectors(N,result_coo,result_crs))
    {
        std::cout<<"The results are not the same!"<<std::endl;
        //std::cout<<"Vectors: "<<std::endl;
        //PrintVector(N,result_coo);
        //PrintVector(N,result_crs);
    }

    delete[] result_coo;
    delete[] result_crs;
    delete[] vector;
    
    return(0);
}
