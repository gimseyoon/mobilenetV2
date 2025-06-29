#include <stdio.h>
#include <stdlib.h>

#define ZP 1
#define S 1
#define X_ROW 512
#define X_COLUMN 512
#define W_ROW 3
#define W_COLUMN 3
#define Y_ROW ((X_ROW -W_ROW + 2*ZP)/S +1)
#define Y_COLUMN ((X_COLUMN -W_COLUMN + 2*ZP)/S +1)
#define THRESHOLD 50

int main()
{   
    ///////////////////////////////////////////////////
    // declare ////////////////////////////////////////
    ///////////////////////////////////////////////////      
    unsigned char x[X_ROW][X_COLUMN]={0,};
    unsigned char x_padded[X_ROW+2*ZP][X_COLUMN+2*ZP] = {0, };
    int w[W_ROW][W_COLUMN] = {{1,0,-1},{1,0,-1},{1,0,-1}};
    int y[Y_ROW][Y_COLUMN] = {0, };
    unsigned char y_after_threshold[Y_ROW][Y_COLUMN] = {0,};



    ///////////////////////////////////////////////////
    // read file (fp -> x ) ///////////////////////////
    ///////////////////////////////////////////////////        
    FILE *fp_input = fopen("C:/seyoon/visual_studio_code/input_output_file/input/lenna8bit.raw","rb");
    if(fp_input == NULL){
        printf("file open error!");
    } 
    fread ((void*)x, sizeof(unsigned char), X_ROW*X_COLUMN, fp_input );
    fclose(fp_input);

    /* print x */
    /*
    for(int i=0; i < X_ROW; i++) 
    {
        for(int j = 0; j < X_COLUMN; j++) 
        {
            printf("x[%d][%d] = %d\n", i,j,x[i][j]);

            if(j==(X_COLUMN - 1)) printf("\n");
        }
    }
    */



    ///////////////////////////////////////////////////
    // padding ( x -> x_padded) ///////////////////////
    ///////////////////////////////////////////////////     
    
    printf("[ padding part ( x -> x_padded) ]\n");
    for(int i=0; i < X_ROW; i++) 
    {
        for(int j = 0; j < X_COLUMN; j++) 
        {
           x_padded[i+1][j+1] = x[i][j]; 
        }
    }
    /* print x_padded */
    /*
    for(int i=0; i < X_ROW + 2; i++) 
    {
        for(int j = 0; j < X_COLUMN + 2; j++) 
        {
            printf("x_padded[%d][%d] = %d\n", i,j,x_padded[i][j]);

            if(j==X_ROW + 1) printf("\n");
        }
    }
    */


    ///////////////////////////////////////////////////
    // compute ////////////////////////////////////////
    ///////////////////////////////////////////////////     

    printf("[ compute part : x_padded ]\n");
    for(int i=0; i < Y_ROW; i++) 
    {
        for(int j = 0; j < Y_COLUMN; j++)
        {
            for(int k = 0; k < W_ROW; k++) 
            {
                for(int p = 0; p < W_COLUMN; p++) 
                {
                    y[i][j] += (x_padded[i+k][j+p]) * (w[k][p]);
                    /*
                    printf("i=%d, j=%d, k=%d, p=%d\t",i,j,k,p);
                    printf("y[%d][%d] = %d\n", i,j,y[i][j]);
                    */
                  
                }     
            }  
        }
    }
    
    printf("\n");



    ///////////////////////////////////////////////////
    // threshold //////////////////////////////////////
    /////////////////////////////////////////////////// 
    for(int i=0; i < Y_ROW; i++) // 0~9
    {
        for(int j = 0; j < Y_COLUMN; j++) // 0~9
        {
            if( abs(y[i][j]) < THRESHOLD) y_after_threshold[i][j] = 0;
            else y_after_threshold[i][j] = 255;
        }
    }


    ///////////////////////////////////////////////////
    // write file (y -> fp ) //////////////////////////
    ///////////////////////////////////////////////////        
    FILE *fp_out = fopen("C:/seyoon/visual_studio_code/input_output_file/output/lenna8bit_output.raw","wb");
    if(fp_out == NULL){
        printf("file open error!");
    } 

    fwrite ((void*)y_after_threshold, sizeof(unsigned char), Y_ROW*Y_COLUMN, fp_out );
    fclose(fp_out);



    ///////////////////////////////////////////////////
    //print result ////////////////////////////////////
    /////////////////////////////////////////////////// 
    
    printf("[ print output result : x_padded ]\n");
    for(int i=0; i < Y_ROW; i++) // 0~9
    {
        for(int j = 0; j < Y_COLUMN; j++) // 0~9
        {
            printf("y[%d][%d] = %d\n", i,j,y_after_threshold[i][j]);
        }
    }


    return 0;
}