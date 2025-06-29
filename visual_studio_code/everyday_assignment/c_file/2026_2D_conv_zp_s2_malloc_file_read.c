#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define X_ROW 10
#define X_COLUMN 10
#define W_ROW 3
#define W_COLUMN 3




int main()
{



    ///////////////////////////////////////////////////
    //declare//////////////////////////////////////////
    /////////////////////////////////////////////////// 
    int x[10][10] = {0, };
    

    //malloc
    int *px = (int *)malloc(sizeof(int)* X_ROW * X_COLUMN );
    if(px==NULL){
        printf("\n[ Memory Allocation error! ]");
        return 1;
    }
    else{
        printf("\n[ Memory Allocation success! ]\n");
    }
    
    //memcpy
    memcpy(px, x, sizeof(x));




    int w[3][3] = {{1,0,-1},{1,0,-1},{1,0,-1}};
    int y[5][5] = {0, };

    int x_padded[12][12] = {0, };
    int x_padded_row = 12, x_padded_column = 12;

    int *pw = &w[0][0];
    int *py = &y[0][0];
    int *px_padded = &x_padded[0][0];


    FILE *fp = fopen("C:/seyoon/visual_studio_code/test_input.txt","r");
    if(fp==NULL){
        printf("error");
        return 1;
    }

    // fscanf로 px에 100개 값 저장
    for (int i = 0; i < X_ROW * X_COLUMN; i++) {
        fscanf(fp, "%d", &px[i]);
    }
    fclose(fp);



    ///////////////////////////////////////////////////
    //padding part ( x -> x_padded)////////////////////
    ///////////////////////////////////////////////////     
    
    printf("\n[ padding start... ]\n");
    for(int i=0; i < ( x_padded_row - W_ROW + 1); i++) // 0~9
    {
        for(int j = 0; j < ( x_padded_column - W_ROW + 1); j++) // 0~9
        {
           *(px_padded + (x_padded_row*(i+1)+(j+1)) ) = *(px + X_ROW*i+j); 
        }
    }
    printf("[ padding complete ]\n\n");








    ///////////////////////////////////////////////////
    //print x_padded///////////////////////////////////
    ///////////////////////////////////////////////////   

    printf("[ printf x_padded ]\n");
    for(int i=0; i < ( x_padded_row); i++) // 0~11
    {
        for(int j = 0; j < ( x_padded_column); j++) // 0~11
        {
            printf("x_p[%d][%d] = %d ", i,j,*(px_padded + ((x_padded_row)*(i)+(j)) ) );

            if(j==(x_padded_row - 1)) printf("\n");
        }
    }
    printf("\n");









    ///////////////////////////////////////////////////
    // compute part : 2D_convolution //////////////////
    ///////////////////////////////////////////////////     

    printf("[ compute part : 2D convolution... ]\n");
    for(int i=0; i < ( (x_padded_row - W_ROW)/2 + 1 ); i++) // 0~4
    {
        for(int j = 0; j < ( (x_padded_column - W_COLUMN)/2  + 1 ); j++) // 0~4
        {
            for(int k = 0; k < W_ROW; k++) // 0~2
            {
                for(int p = 0; p < W_COLUMN; p++) // 0~2
                {
                    *(py + ((x_padded_column - W_COLUMN)/2  + 1 )*i+j) += (*(px_padded + ((x_padded_row)*(2*i+k)+(2*j+p)) )) * (*(pw + W_ROW*k+p));
                    
                    printf("i=%d, j=%d, k=%d, p=%d\t",i,j,k,p);
                    printf("y[%d][%d] = %d\n", i,j,*(py + ((x_padded_column - W_COLUMN)/2  + 1 )*i+j));
                    
                  
                }     
            }  
        }
    }
    printf("[ compute complete ]\n\n");







    ///////////////////////////////////////////////////
    //print result : x_padded /////////////////////////
    /////////////////////////////////////////////////// 
    
    printf("[ print result ]\n");

    
    for(int i=0; i < ( (x_padded_row - W_ROW)/2 + 1 ); i++) // 0~4
    {
        for(int j = 0; j < ( (x_padded_column - W_ROW)/2 + 1); j++) // 0~4
        {
            printf("y[%d][%d] = %d\t", i,j,*(py + ((x_padded_column - W_COLUMN)/2  + 1 )*i+j));

            if(j==(x_padded_row - W_ROW -1)/2 ) printf("\n");
        }
    }

    free(px);

    return 0;
} 
