/*
 full connection 
 - input: 10
 - output: 5
 - weight: 5x10

*/


#include <stdio.h>

#define INPUT_SIZE 10
#define W_ROW 10
#define W_COLUMN 5
#define OUTPUT_SIZE 5

int main()
{

    ///////////////////////////////////////////////////
    // declare ////////////////////////////////////////
    /////////////////////////////////////////////////// 
    int x[INPUT_SIZE] = {0,};
    int w[W_ROW][W_COLUMN] = {0, };
    int y[OUTPUT_SIZE] = {0,};


    ///////////////////////////////////////////////////
    // input file READ ////////////////////////////////
    ///////////////////////////////////////////////////     

    FILE *fp_i = fopen("C:/seyoon/visual_studio_code/everyday_assignment/i_o_w_file/input/full_connection_input.txt", "rb");
    if (fp_i == NULL) {
        printf("file read error!\n");
        return 1;
    }    
    for (int i = 0; i < INPUT_SIZE; i++) {
        fscanf(fp_i, "%d", &x[i]);
    }
    fclose(fp_i);
    
    for (int i = 0; i < INPUT_SIZE; i++) {
        printf("x[%d] = %d\n", i, x[i]);
    }
    

    ///////////////////////////////////////////////////
    // weight file READ ///////////////////////////////
    ///////////////////////////////////////////////////  

    FILE *fp_w = fopen("C:/seyoon/visual_studio_code/everyday_assignment/i_o_w_file/weight/full_connection_w.txt", "rb");
    if (fp_w == NULL) {
        printf("file read error!\n");
        return 1;
    }    
    for (int i = 0; i < W_ROW; i++) {
        for(int j = 0; j < W_COLUMN; j++){
            fscanf(fp_w, "%d", &w[i][j]);
        }  
    }
    fclose(fp_w);
    
    for (int i = 0; i < W_ROW; i++) {
        for(int j = 0; j < W_COLUMN; j++){
            printf("w[%d][%d] = %d\t", i,j,w[i][j]);
            if(j==W_COLUMN-1) printf("\n");
        }  
    }
    
    ///////////////////////////////////////////////////
    // compute ////////////////////////////////////////
    /////////////////////////////////////////////////// 

    for(int i = 0; i<W_COLUMN; i++)
    {
        for(int j = 0; j<W_ROW; j++)
        {
            y[i] += x[j]*w[j][i];
        }
    }

    ///////////////////////////////////////////////////
    // print result ///////////////////////////////////
    /////////////////////////////////////////////////// 
    for (int i = 0; i < OUTPUT_SIZE; i++) {
        printf("y[%d] = %d\n", i, y[i]);
    }    


    return 0;
}
   