/***************************************************************************************
 * 
 * t3.c: producto tensorial calculado con mpi y thread
 *
 * Programmer: Cristobal Gallardo Cubillos & Vicente Santos Varas
 *
 * Santiago de Chile, 7/12/2023
 *
 **************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include "mpi.h"
#include <pthread.h>

#define MASTER 0
#define blockxthread 16
MPI_Status status;
float **Matrix1, **Matrix2, **Resultado, **Matrix1t, **Matrix2t;
unsigned int num_threads, partition_type, rt, ct, rbt, cbt;
/*
 *
 */
struct data {
    unsigned int id; 
};

/*
 *
 */
void Usage(char *mess) {

    printf("\nUsage: mpirun -np k %s -P -O data.txt\n",mess);
    printf("K = numero de procesos\n");
    printf("P = {V: particion vertical, H: particion horizontal}\n");
    printf("O = {S: modo silencioso, V: modo vervoso}\n\n");
}

/*
 *
 */ //llena la matrix
float **FillMatrix(unsigned int r, unsigned int c, FILE *archivo) {
    unsigned int i, j;
    float **mat;

    mat = calloc(r, sizeof(float *));
    for (i = 0; i < r; i = i + 1) {
        mat[i] = calloc(c, sizeof(float));
    }
    for (i = 0; i < r; i = i + 1) {
        for (j = 0; j < c; j = j + 1) {
            fscanf(archivo, "%f", &mat[i][j]);
        }
    }
    return mat;
}

/*
 *
 */
void PrintMatrix(unsigned int r, unsigned int c, float **mat) {
    
   unsigned int i, j;
   
   for (i = 0; i < r; i = i + 1) {
      for (j = 0; j < c; j = j + 1)
         printf(" %.2f ",mat[i][j]);
      printf("\n");
   }
}

/*
 *
 */ //proceso de cálculo
float **Process(float **x, float **y, unsigned int r, unsigned c,unsigned rb, unsigned int cb) {

   unsigned int i, j, startRow, startCol, k, l;
   float **res;

   res = calloc(r * rb,sizeof(float* ));
   for (i = 0; i < r * rb; i = i + 1)
      res[i] = calloc(c * cb,sizeof(float*));
   for(i = 0; i < r; i = i + 1){
        for(j = 0; j < c; j = j + 1){
            startRow = i * rb;
            startCol = j * cb;
            for(k = 0; k < rb; k = k + 1){
                for(l = 0; l < cb; l = l + 1){
                    res[startRow + k][startCol + l] = x[i][j] * y[k][l];
                }
            }
        }
    }
    return res;
}

/*
 *
 */
 void printMatrix1D(float *matrix, int m, int k){
    printf("Matriz aplanada:\n");
    int i;
    for ( i = 0; i < m * k; i = i + 1){
        printf("%.2f ", matrix[i]);
    }
    printf("\n");
}

/*
 *
 */
void FreeMatrix(float **mat, unsigned int r) {

   unsigned int i;
   
   for (i = 0; i < r; i = i + 1)
      free(mat[i]);
   free(mat);
}

/*
 *
 */
float **GetSubMatrix(float **y, unsigned int start_row, unsigned int start_col, unsigned int end_row, unsigned int end_col) {
    float** x;
    unsigned int i, j;
    x = (float **)malloc((end_row - start_row ) * sizeof(float *));
    for (i = 0; i < (end_row - start_row); i = i + 1) {
        x[i] = (float *)malloc((end_col - start_col ) * sizeof(float));
    }
    for (i = start_row; i < end_row; i = i + 1) {
        for (j = start_col; j < end_col; j = j + 1) {
            x[i - start_row][j - start_col] = y[i][j];
        }
    }
    return x;
}

/*
 *
 */
float* Matrix1D(float **matrix, int filas, int columnas){
    float *matriz_aplanada = (float*) malloc(filas * columnas * sizeof(float));
    int i, j, indice = 0;

    for (i = 0; i < filas; i = i + 1){    
        for ( j = 0; j < columnas; j = j + 1){
            matriz_aplanada[indice] =matrix[i][j];
            indice = indice + 1;
        }
    }
    return matriz_aplanada;
}

/*
 *
 */
float** Matrix2D(float* matriz_aplanada , int filas, int columnas){
    float** matrix = (float**) malloc(filas * sizeof(float*));
    for (int i = 0; i < filas; i = i + 1) {
    matrix[i] = (float*) malloc(columnas * sizeof(float));
    }

    int indice = 0, i, j;
    for (i = 0; i < filas; i = i + 1){
        for ( j = 0; j < columnas; j = j + 1){
            matrix[i][j]=matriz_aplanada[indice];
            indice = indice + 1;
        }
    }
    return matrix;
}

/*
 *
 */
 void* tensorial_product(void* args) {
    struct data *data = (struct data*)args;
    int id = data->id;
    unsigned int i, j, startRow, startCol, k, l;
    
   
    if ( partition_type == 0) { // Partición horizontal
        for(i = id; i < rt; i = i + num_threads){
            for(j = 0; j < ct; j = j + 1){
                startRow = i * rbt;
                startCol = j * cbt;
                for(k = 0; k < rbt; k = k + 1){
                    for(l = 0; l < cbt; l = l + 1){
                        Resultado[startRow + k][startCol + l] = Matrix1t[i][j] * Matrix2t[k][l];
                    }
                }
            }
        }
        
    } else if ( partition_type == 1) { // Partición vertical
        for(j = id; j < ct; j = j + num_threads){
            for(i = 0; i < rt; i = i + 1){
                startRow = i * rbt;
                startCol = j * cbt;
                for(k = 0; k < rbt; k = k + 1){
                    for(l = 0; l < cbt; l = l + 1){
                        Resultado[startRow + k][startCol + l] = Matrix1t[i][j] * Matrix2t[k][l];
                    }
                }
            }
        }
    }
    
    pthread_exit(NULL);
}

/*
 *
 */
int main(int argc, char **argv) {
    int m, i, k, n, ii, start_row, end_row, start_col, end_col, p = 0, o = 0, rows_per_process, colum_per_process, extra_rows, extra_colum, rank, size, me, datos[8], rows_per_processb, extra_rowsb, end_rowb, start_rowb, colum_per_processb, extra_columb, start_colb, end_colb, sentmess = 0, recvmess = 0, j; 
    float* matrixA_1D, *matrixB_1D, *matrixC_1D, **matrixC, **arreglo, **arreglob, E_cpu;
    char  processor_name[MPI_MAX_PROCESSOR_NAME], nombre[20]; 
    long E_wall;  
    FILE *archivo;
    time_t  ts, te;
    clock_t cs, ce;


    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(processor_name, &me);

    if(rank == MASTER){   
        ts = time(NULL);
        cs = clock();
        if (strcmp(argv[1],"-V") == 0) {       
            strcpy(nombre, "Vertical");
        }
        else if (strcmp(argv[1],"-H") == 0){
            strcpy(nombre, "Horizontal");
        }
        if (argc == 4) {
            // Leer las dimensiones de la matriz
            archivo = fopen(argv[3],"r");
            if (archivo == NULL) {
                printf("No se pudo abrir el archivo.\n");
                return 0;
            }
            fscanf(archivo, "%d", &m); // Filas A
            fscanf(archivo, "%d", &k); // Columnas A y filas B
            fscanf(archivo, "%d", &n); // Columnas B
            datos[1] = k;
            datos[2] = n;
            
            printf("Numero de Procesos: %d, Modo particion: %s, Dimension Matriz A[%d,%d] y Matriz B[%d,%d]\n", size, nombre, m, k, k, n);  
            fflush(stdout);
            // Leer e inicializa las matrices
            Matrix1 = FillMatrix(m, k, archivo);
            Matrix2 = FillMatrix(k, n, archivo);
            fclose(archivo);
            
            
            
            // Muestra las matrices si se desean ver
            if (strcmp(argv[2], "-V") == 0) {
                datos[5] = 1;
                printf(" Matriz A(%d,%d):\n\n", m, k);
                PrintMatrix(m,k,Matrix1);
                printf("\n");
                printf(" Matriz B(%d,%d):\n\n", k, n);
                PrintMatrix(k, n, Matrix2);
            }
            //Pasar la matrix de 2D a 1D
            matrixA_1D = Matrix1D(Matrix1, m, k);
            matrixB_1D = Matrix1D(Matrix2, k, n);

            if (strcmp(argv[1], "-H") == 0) {       //Modo horizontal
                rows_per_process = m / (size);
                rows_per_processb = k / (size);
                extra_rows = m % (size);
                extra_rowsb = k % (size);
                datos[0] = rows_per_process;
                datos[3] = k;
                partition_type = 0;
                datos[7] = 0;
                
                if (strcmp(argv[2], "-V") == 0) //muestra las Submatrices
                {
                    for (i = 0; i < size; i = i + 1){ 
                        start_row = p;
                        start_rowb = o;
                        end_row = p + rows_per_process;
                        end_rowb = o + rows_per_processb;

                        if(extra_rows > 0 && i == size - 1){             // Caso que no quepan todas las filas de a en los nodos asignados las que sobren se llevaran al ultimo nodo
                            end_row = end_row + extra_rows;
                            rows_per_process = rows_per_process + extra_rows;
                            
                        }
                        if(extra_rowsb > 0 && i == size - 1){             // Caso que no quepan todas las filas de b en los nodos asignados las que sobren se llevaran al ultimo nodo
                            end_rowb = end_rowb + extra_rowsb;
                            rows_per_processb = rows_per_processb + extra_rowsb;
                            
                        }

                        arreglo = GetSubMatrix(Matrix1, start_row, 0, end_row, k);
                        arreglob = GetSubMatrix(Matrix2, start_rowb, 0, end_rowb, n);
                        if(strcmp(argv[2], "-V") == 0){
                            if (i == 0){
                                printf("\nSubMatrix A en el proceso MASTER\n");      
                            }
                            else{
                                printf("\nSubMatrix A en el proceso %d\n", i);
                            }
                            PrintMatrix(rows_per_process, k, arreglo); // Imprime las submatrices de A y B enviadas a cada proceso
                            fflush(stdout);
                            if (i == 0){
                                printf("\nSubMatrix B en el proceso MASTER\n");      
                            }
                            else{
                                printf("\nSubMatrix B en el proceso %d\n", i);
                            }
                            PrintMatrix(rows_per_processb, n, arreglob);
                            fflush(stdout);
                        }
                        p = p + rows_per_process;
                        o = o + rows_per_processb;
                    }
                }
                if(strcmp(argv[2], "-V") == 0){
                    printf("\nMatriz Resultado\n");
                }
                
                rows_per_process = m / (size);
                rows_per_processb = k / (size);
                extra_rows = m % (size);
                extra_rowsb = k % (size);
                datos[0] = rows_per_process;
                datos[3] = k;
                datos[1] = rows_per_processb;
                datos[2] = n;
                datos[4] = rows_per_processb + extra_rowsb;
                datos[6] = n;
                p = 0;
                o = 0;
                num_threads = rows_per_process;
                pthread_t threads[num_threads];
                struct data data[num_threads];
                rt = rows_per_process;
                ct = k;
                rbt = rows_per_processb;
                cbt = n;
                
                for (i = 0; i < size; i = i + 1){ 
                    start_row = p;
                    start_rowb = o;
                    end_row = p + rows_per_process;
                    end_rowb = o + rows_per_processb;

                    if(extra_rows > 0 && i == size - 1){             //Caso que no quepan todas las filas de a en los nodos asignados las que sobren se llevaran al ultimo nodo
                        end_row = end_row + extra_rows;
                        rows_per_process = rows_per_process + extra_rows;
                        datos[0] = rows_per_process;
                    }
                        
                    if(extra_rowsb > 0 && i == size - 1){             //Caso que no quepan todas las filas de b en los nodos asignados las que sobren se llevaran al ultimo nodo
                        end_rowb = end_rowb + extra_rowsb;
                        rows_per_processb = rows_per_processb + extra_rowsb;
                        datos[1] = rows_per_processb;
                    }
                    arreglo = GetSubMatrix(Matrix1, start_row, 0, end_row, k);
                    arreglob = GetSubMatrix(Matrix2, start_rowb, 0, end_rowb, n);
                    Resultado = (float **)calloc(rows_per_process * rows_per_processb, sizeof(float *));
                    for (ii = 0; ii < rows_per_process * rows_per_processb; ii = ii + 1) {
                        Resultado[ii] = (float *)calloc(k * n, sizeof(float));
                    }
                    if (i == 0){
                       if (rows_per_process > 1){
                            Matrix1t = arreglo;
                            Matrix2t = arreglob;
                            for ( j = 0; j < num_threads; j = j + 1) {      //Se crean hilos en funcion al numero de hebras predefinido por la cantidad de filas 
                                data[j].id = j;
                                pthread_create(&threads[j], NULL, tensorial_product, &data[j]);
                            }
                            for (j = 0; j < num_threads; j = j + 1) {       //Se asegura que terminen todos los hilos su ejecucion.
                                pthread_join(threads[j], NULL);
                            }
                        }
                        else{
                            Resultado = Process(arreglo, arreglob, rows_per_process, k, rows_per_processb, n);
                        }
                        
                        
                        if(strcmp(argv[2], "-V") == 0){
                            PrintMatrix(rows_per_process * rows_per_processb, k * n, Resultado);
                        }                
                    }
                    else{
                        matrixA_1D = Matrix1D(arreglo, rows_per_process, k);
                        matrixB_1D = Matrix1D(arreglob, rows_per_processb,n);
                        // Mandar matrices a los nodos
                        MPI_Send(datos , 8, MPI_INT, i, i, MPI_COMM_WORLD);
                        MPI_Send(matrixA_1D, rows_per_process * k, MPI_FLOAT, i, i, MPI_COMM_WORLD);   //Se envian los datos y las matrices
                        MPI_Send(matrixB_1D , rows_per_processb * n, MPI_FLOAT, i, i , MPI_COMM_WORLD);
                        sentmess = sentmess + 3;
                    }
                    p = p + rows_per_process;
                    o = o + rows_per_processb;
                } 
                
                rows_per_process = m / (size);
                rows_per_processb = k / (size);
                extra_rows = m % (size);
                extra_rowsb = k % (size);
                p = 0;
                o = 0;
                datos[1] = rows_per_processb;
                arreglo = GetSubMatrix(Matrix1, p, 0, rows_per_process, k);
                arreglob = GetSubMatrix(Matrix2, o, 0, rows_per_processb, n);
                matrixB_1D = Matrix1D(arreglob, rows_per_processb,n);
                        
                for (i = 1; i < size; i = i + 1){   //Proceso que recibe los datos de matriz B desde procesos y calcula la matriz restante para el proceso MASTER
                    start_row = p;
                    start_rowb = o;
                    end_row = p + rows_per_process;
                    end_rowb = o + rows_per_processb;
                    
                    MPI_Send(datos, 8, MPI_INT, i, 0, MPI_COMM_WORLD);
                    MPI_Send(matrixB_1D, rows_per_processb * n, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
                    sentmess = sentmess + 2;
                    
                    if(extra_rowsb > 0 && i == size - 1){             //Caso que no quepan todas las filas de b en los nodos asignados las que sobren se llevaran al ultimo nodo
                        end_rowb = end_rowb + extra_rowsb;
                        rows_per_processb = rows_per_processb + extra_rowsb;
                    }
                    matrixC_1D = (float*) malloc(rows_per_processb * n * sizeof(float));
                    MPI_Recv(datos , 8, MPI_INT, i, i, MPI_COMM_WORLD, &status); 
                    MPI_Recv(matrixC_1D, rows_per_processb * n, MPI_FLOAT, i, i, MPI_COMM_WORLD, &status);
                    recvmess = recvmess + 2;
                    rt = rows_per_process;
                    rbt = rows_per_processb;
                
                    
                    matrixC = Matrix2D(matrixC_1D, rows_per_processb, n);
                    if (rows_per_process > 1){
                        Matrix1t = arreglo;
                        Matrix2t = matrixC;
                        for ( j = 0; j < num_threads; j = j + 1) {
                            data[j].id = j;
                            pthread_create(&threads[j], NULL, tensorial_product, &data[j]);
                        }
                        for (j = 0; j < num_threads; j = j + 1) {       //Se asegura que terminen todos los hilos su ejecucion.
                            pthread_join(threads[j], NULL);
                        }
                    }
                    else{
                        Resultado = Process(arreglo, matrixC, rows_per_process, k, rows_per_processb, n);
                    }
                    
                    if(strcmp(argv[2], "-V") == 0){
                        PrintMatrix(rows_per_process * rows_per_processb, k * n, Resultado);
                    }
                    
                    o = o + rows_per_processb;
                }
            }
            else{ //Modo vertical
                partition_type = 1;
                datos[7] = 1;
                colum_per_process = k / (size);
                colum_per_processb = n / (size);
                extra_colum = k % (size);
                extra_columb = n % (size);
                datos[0] = m;
                datos[3] = colum_per_process;
                
                if (strcmp(argv[2], "-V") == 0) //muestra las Submatrices
                {
                    for (i = 0; i < size; i = i + 1){ 
                        start_col = p;
                        start_colb = o;
                        end_col = p + colum_per_process;
                        end_colb = o + colum_per_processb;
                        if(extra_colum > 0 && i == size - 1){             //Caso que no quepan todas las columnas de A en los nodos asignados las que sobren se llevaran al ultimo nodo
                            end_col = end_col + extra_colum;
                            colum_per_process = colum_per_process + extra_colum;
                        }
                        
                        if(extra_columb > 0 && i == size - 1){             //Caso que no quepan todas las columnas de B en los nodos asignados las que sobren se llevaran al ultimo nodo
                            end_colb = end_colb + extra_columb;
                            colum_per_processb = colum_per_processb + extra_columb;
                        }
                        
                        arreglo = GetSubMatrix(Matrix1, 0, start_col, m, end_col);
                        arreglob = GetSubMatrix(Matrix2, 0, start_colb, k, end_colb);
                        if(strcmp(argv[2], "-V") == 0){
                            if (i == 0){
                                printf("\nSubMatrix A en el proceso MASTER\n");      
                            }
                            else{
                                printf("\nSubMatrix A en el proceso %d\n", i);
                            }
                            PrintMatrix(m, colum_per_process, arreglo);     //Imprime las submatrices enviadas a cada proceso
                            fflush(stdout);
                            if (i == 0){
                                printf("\nSubMatrix B en el proceso MASTER\n");      
                            }
                            else{
                                printf("\nSubMatrix B en el proceso %d\n", i);
                            }
                            PrintMatrix(k, colum_per_processb, arreglob);
                            fflush(stdout);
                        }
                        p = p + colum_per_process;
                        o = o + colum_per_processb;
                    }
                }
                if(strcmp(argv[2], "-V") == 0){
                    printf("\nMatriz Resultado\n");
                }
                colum_per_process = k / (size);
                colum_per_processb = n / (size);
                extra_colum = k % (size);
                extra_columb = n % (size);
                num_threads = colum_per_process;
                pthread_t threads[num_threads];
                struct data data[num_threads];
                p = 0;
                o = 0;
                datos[0] = m;
                datos[3] = colum_per_process;
                datos[1] = k;
                datos[2] = colum_per_processb;
                datos[4] = k;
                datos[6] = colum_per_processb + extra_columb;
                rt = m;
                ct = colum_per_process;
                rbt = k;
                cbt = colum_per_processb;

                for (i = 0; i < size; i = i + 1){ 
                    start_col = p;
                    start_colb = o;
                    end_col = p + colum_per_process;
                    end_colb = o + colum_per_processb;

                    if(extra_colum > 0 && i == size - 1){             //Caso que no quepan todas las columnas de a en los nodos asignados las que sobren se llevaran al ultimo nodo
                        end_col = end_col + extra_colum;
                        colum_per_process = colum_per_process + extra_colum;
                        datos[3] = colum_per_process;
                    }
                        
                    if(extra_columb > 0 && i == size - 1){             //Caso que no quepan todas las columnas de b en los nodos asignados las que sobren se llevaran al ultimo nodo
                        end_colb = end_colb + extra_columb;
                        colum_per_processb = colum_per_processb + extra_columb;
                        datos[2] = colum_per_processb;
                    }
                    arreglo = GetSubMatrix(Matrix1, 0, start_col, m, end_col);
                    arreglob = GetSubMatrix(Matrix2, 0, start_colb, k, end_colb);
                    Resultado = (float **)calloc(m * k, sizeof(float *));
                        for (ii = 0; ii < m * k; ii = ii + 1) {
                             Resultado[ii] = (float *)calloc(colum_per_process * colum_per_processb, sizeof(float));
                        }
                    if (i == 0){
                        if (colum_per_process > 1){
                            Matrix1t = arreglo;
                            Matrix2t = arreglob;
                            for ( j = 0; j < num_threads; j = j + 1) {      //Se crean hilos en funcion a la cantidad de columnas
                                data[j].id = j;
                                pthread_create(&threads[j], NULL, tensorial_product, &data[j]);
                            }
                            for (j = 0; j < num_threads; j = j + 1) {       //Se asegura que terminen todos los hilos su ejecucion.
                                pthread_join(threads[j], NULL);
                            }
                        }
                        else{
                            Resultado = Process(arreglo, arreglob, m, colum_per_process, k, colum_per_processb);
                        }
                        if(strcmp(argv[2], "-V") == 0){
                            PrintMatrix(m * k, colum_per_process * colum_per_processb, Resultado);
                        }              
                    }
                    else{
                        matrixA_1D = Matrix1D(arreglo, m, colum_per_process);
                        matrixB_1D = Matrix1D(arreglob, k, colum_per_processb);
                        // Mandar matrices a los nodos
                        MPI_Send(datos , 8, MPI_INT, i, i, MPI_COMM_WORLD);
                        MPI_Send(matrixA_1D, m * colum_per_process, MPI_FLOAT, i, i, MPI_COMM_WORLD);   //Se envian los datos y las matrices
                        MPI_Send(matrixB_1D , k * colum_per_processb, MPI_FLOAT, i, i , MPI_COMM_WORLD);
                        sentmess = sentmess + 3;
                    }
                    p = p + colum_per_process;
                    o = o + colum_per_processb;
                } 
                
                colum_per_process = k / (size);
                colum_per_processb = n / (size);
                extra_colum = k % (size);
                extra_columb = n % (size);
                p = 0;
                o = 0;
                datos[2] = colum_per_processb;
                arreglo = GetSubMatrix(Matrix1, 0, p, m, end_col);
                arreglob = GetSubMatrix(Matrix2, 0, o, k, end_colb);
                matrixB_1D = Matrix1D(arreglob, k, colum_per_processb);
                        
                for (i = 1; i < size; i = i + 1){   //Proceso que recibe los datos de matriz B desde procesos y calcula la matriz restante para el proceso MASTER
                    start_col = p;
                    start_colb = o;
                    end_col = p + colum_per_process;
                    end_colb = o + colum_per_processb;
                    
                    MPI_Send(datos, 8, MPI_INT, i, 0, MPI_COMM_WORLD);
                    MPI_Send(matrixB_1D, k * colum_per_processb, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
                    sentmess = sentmess + 2;
                    
                    if(extra_columb > 0 && i == size - 1){             //Caso que no quepan todas las filas en los nodos asignados las que sobren se llevaran al ultimo nodo
                        end_colb = end_colb + extra_columb;
                        colum_per_processb = colum_per_processb + extra_columb;
                    }
                    matrixC_1D = (float*) malloc(k * colum_per_processb * sizeof(float));
                    MPI_Recv(datos , 8, MPI_INT, i, i, MPI_COMM_WORLD, &status); 
                    MPI_Recv(matrixC_1D, k * colum_per_processb, MPI_FLOAT, i, i, MPI_COMM_WORLD, &status);
                    recvmess = recvmess + 2;
                    ct = colum_per_process;
                    cbt = colum_per_processb;
                    
                    matrixC = Matrix2D(matrixC_1D, k, colum_per_processb);
                    if (colum_per_process > 1){
                        Matrix1t = arreglo;
                        Matrix2t = matrixC;
                        for ( j = 0; j < num_threads; j = j + 1) {
                            data[j].id = j;
                            pthread_create(&threads[j], NULL, tensorial_product, &data[j]);
                        }
                        for (j = 0; j < num_threads; j = j + 1) {       //Se asegura que terminen todos los hilos su ejecucion.
                            pthread_join(threads[j], NULL);
                        }
                    }
                    else{
                        Resultado = Process(arreglo, matrixC, m, colum_per_process, k, colum_per_processb);
                    }
                    
                    if(strcmp(argv[2], "-V") == 0){
                        PrintMatrix(m * k, colum_per_process * colum_per_processb, Resultado);
                    }
                    
                    o = o + colum_per_processb;
                }
                
            }
            ce = clock();
            te = time(NULL);
            E_wall = (long) (te - ts);
            E_cpu = (float)(ce - cs) / CLOCKS_PER_SEC;
            MPI_Barrier( MPI_COMM_WORLD);
            printf("From MASTER - Elapsed CPU Time %f Wall Time %ld \n\t Sent Mess: %d, Recv Mess: %d \n", E_cpu, E_wall, sentmess, recvmess); 
        }
        else {
            Usage(argv[0]);
        }
    }
    else{//Procesos 
        ts = time(NULL); //inicia el contador de tiempo=0;
        cs = clock();
        
        MPI_Recv(datos , 8, MPI_INT, MASTER , rank , MPI_COMM_WORLD, &status); //Se reciben los datos
        rows_per_process = datos[0];
        rows_per_processb = datos[1];
        colum_per_processb = datos[2];
        colum_per_process = datos[3];
        partition_type = datos[7];

        num_threads = rows_per_process;
        pthread_t threads[num_threads];
        struct data data[num_threads];

        matrixA_1D = (float*) malloc(rows_per_process * colum_per_process * sizeof(float));    //Se le asigna tamaño a las matrices unidimensionales
        matrixB_1D = (float*) malloc(rows_per_processb * colum_per_processb * sizeof(float));
        matrixC_1D = (float*) malloc(rows_per_processb * colum_per_processb * sizeof(float));
        MPI_Recv(matrixA_1D, rows_per_process*colum_per_process, MPI_FLOAT, MASTER, rank, MPI_COMM_WORLD, &status);  //Se reciben las matrices
        MPI_Recv(matrixB_1D, rows_per_processb * colum_per_processb , MPI_FLOAT, MASTER, rank, MPI_COMM_WORLD, &status);
        recvmess = recvmess + 3;

        Matrix1 = Matrix2D(matrixA_1D, rows_per_process, colum_per_process);
        Matrix2 = Matrix2D(matrixB_1D, rows_per_processb, colum_per_processb);

        
        for  (i = 0; i < size; i = i + 1){
            if(i != rank){
                MPI_Send(datos, 8, MPI_INT, i, rank, MPI_COMM_WORLD);
                MPI_Send(matrixB_1D, rows_per_processb * colum_per_processb, MPI_FLOAT, i, rank, MPI_COMM_WORLD); //Se manda la matriz B a los otros procesos
                sentmess = sentmess + 2;
             }   
        }
        for (i = 0; i < size; i = i + 1){
           
            if(i == size - 1){
                rows_per_processb = datos[4];
                colum_per_processb = datos[6];
            }
            Resultado = (float **)calloc(rows_per_process * rows_per_processb, sizeof(float *)); //Se le asigna memoria dinamica a la matriz resultado
            for (ii = 0; ii < rows_per_process * rows_per_processb; ii = ii + 1) {
                Resultado[ii] = (float *)calloc(colum_per_process * colum_per_processb, sizeof(float));
            }
            
            if(i == rank){
                rt = rows_per_process;
                rbt = rows_per_processb;
                ct = colum_per_process;
                cbt = colum_per_processb;
                
                if (rows_per_process > 1 && partition_type == 0 || colum_per_process > 1 && partition_type == 1){ //Se entra a los hilos
                    Matrix1t = Matrix1;
                    Matrix2t = Matrix2;
                    for ( j = 0; j < num_threads; j = j + 1) {
                        data[j].id = j;
                        pthread_create(&threads[j], NULL, tensorial_product, &data[j]);
                    }
                    for (j = 0; j < num_threads; j = j + 1) {       //Se asegura que terminen todos los hilos su ejecucion.
                        pthread_join(threads[j], NULL);
                    }
                        }
                else{
                    Resultado = Process(Matrix1, Matrix2, rows_per_process, colum_per_process, rows_per_processb, colum_per_processb);
                }

                if (datos[5] == 1){   //Si es verboso imprimir la matriz
                    PrintMatrix(rows_per_process * rows_per_processb, colum_per_process * colum_per_processb, Resultado);
                }
            }
            else{    
                MPI_Recv(datos , 8, MPI_INT, i, i , MPI_COMM_WORLD, &status); 
                rows_per_processb = datos[1];
                colum_per_processb = datos[2];
                matrixC_1D = (float*) malloc(rows_per_processb * colum_per_processb * sizeof(float));
                MPI_Recv(matrixC_1D, rows_per_processb * colum_per_processb , MPI_FLOAT, i, i , MPI_COMM_WORLD, &status);
                recvmess = recvmess + 2;
                rt = rows_per_process;
                rbt = rows_per_processb;
                ct = colum_per_process;
                cbt = colum_per_processb;
                
                matrixC = Matrix2D(matrixC_1D, rows_per_processb, colum_per_processb);

                if (rows_per_process > 1 && partition_type == 0 || colum_per_process > 1 && partition_type == 1){
                    Matrix1t = Matrix1; //   Se actualizan los valores de las matrices globales para que puedan trabajar los hilos
                    Matrix2t = matrixC;
                    for ( j = 0; j < num_threads; j = j + 1) {  
                        data[j].id = j;
                        pthread_create(&threads[j], NULL, tensorial_product, &data[j]);
                    }
                    for (j = 0; j < num_threads; j = j + 1) {       //Se asegura que terminen todos los hilos su ejecucion.
                        pthread_join(threads[j], NULL);
                    }
                        }
                else{
                    Resultado = Process(Matrix1, matrixC, rows_per_process, colum_per_process, rows_per_processb, colum_per_processb);
                }

                if (datos[5] == 1){
                    PrintMatrix(rows_per_process * rows_per_processb, colum_per_process * colum_per_processb, Resultado);
                }
                
                
            }
            
        }
       
        ce = clock();
        te = time(NULL);
        E_wall = (long) (te - ts);
        E_cpu = (float)(ce - cs) / CLOCKS_PER_SEC;  //Termina el contador y se muestra tiempo en pantalla.
       
        MPI_Barrier( MPI_COMM_WORLD);               //Barrera para que no se impriman los tiempos en medio de los resultados.
        printf("\nFrom %d - Elapsed CPU Time %f Wall Time %ld \n\t Sent Mess: %d, Recv Mess: %d \n",rank, E_cpu, E_wall, sentmess, recvmess); 
        free(Matrix1);  //Se libera la memoria de las matrices
        free(Matrix2);
        free(matrixC);
        free(matrixA_1D);
        free(matrixB_1D);
        free(matrixC_1D);
        free(Resultado);
    }
    MPI_Finalize();
    return 0;
}