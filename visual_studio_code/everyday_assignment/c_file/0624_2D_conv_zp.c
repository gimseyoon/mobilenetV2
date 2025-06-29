#include <stdio.h>

#define X_ROW 10
#define X_COLUMN 10
#define W_ROW 3
#define W_COLUMN 3




int main()
{


    ///////////////////////////////////////////////////
    //declare//////////////////////////////////////////
    /////////////////////////////////////////////////// 
    int x[10][10] = {{4,5,6,1,2,3,4,1,3,5},
                     {4,5,6,1,2,3,4,1,3,5},
                     {4,5,6,1,2,3,4,1,3,5},
                     {4,5,6,1,2,3,4,1,3,5},
                     {4,5,6,1,2,3,4,1,3,5},
                     {4,5,6,1,2,3,4,1,3,5},
                     {4,5,6,1,2,3,4,1,3,5},
                     {4,5,6,1,2,3,4,1,3,5},
                     {4,5,6,1,2,3,4,1,3,5},
                     {4,5,6,1,2,3,4,1,3,5}};
    
    int w[3][3] = {{1,0,-1},{1,0,-1},{1,0,-1}};
    int y[10][10] = {0, };

    int x_padded[12][12] = {0, };
    int x_padded_row = 12, x_padded_column = 12;




    ///////////////////////////////////////////////////
    //padding part ////////////////////
    ///////////////////////////////////////////////////     
    
    printf("[ padding part ( x -> x_padded) ]\n");
    for(int i=0; i < ( x_padded_row - W_ROW + 1); i++) // need to modify. (i < X_ROW)
    {
        for(int j = 0; j < ( x_padded_column - W_ROW + 1); j++) // need to modify. (j < X_COLUMN)
        {
           x_padded[i+1][j+1] = x[i][j]; 
        }
    }

    printf("\n");


    ///////////////////////////////////////////////////
    // compute part ////////////////////////
    ///////////////////////////////////////////////////     

    printf("[ compute part : x_padded ]\n");
    for(int i=0; i < ( x_padded_row - W_ROW + 1); i++) // 0~9
    {
        for(int j = 0; j < ( x_padded_column - W_COLUMN + 1); j++) // 0~9
        {
            for(int k = 0; k < W_ROW; k++) // 0~2
            {
                for(int p = 0; p < W_COLUMN; p++) // 0~2
                {
                    y[i][j] += (x_padded[i+k][j+p]) * (w[k][p]);
                    
                    printf("i=%d, j=%d, k=%d, p=%d\t",i,j,k,p);
                    printf("y[%d][%d] = %d\n", i,j,y[i][j]);
                    
                  
                }     
            }  
        }
    }
    
    printf("\n");


    ///////////////////////////////////////////////////
    //print result /////////////////////////
    /////////////////////////////////////////////////// 
    
    printf("[ print result : x_padded ]\n\n");
    for(int i=0; i < ( x_padded_row - W_ROW + 1); i++) // 0~9
    {
        for(int j = 0; j < ( x_padded_column - W_ROW + 1); j++) // 0~9
        {
            printf("y[%d][%d] = %d\t", i,j,y[i][j]);

            if(j==(x_padded_row - W_ROW)) printf("\n");
        }
    }



    return 0;
} 
