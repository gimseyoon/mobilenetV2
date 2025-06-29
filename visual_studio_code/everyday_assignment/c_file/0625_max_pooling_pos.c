#include <stdio.h>

#define X_ROW 10         //10 by 10 input
#define POOLING_ROW 2    //2 by 2 pooling

int main(void)
{
    int x[10][10] = {0, };
    int x_pool[5][5] = {0,};
    for(int i = 0; i< 10; i++)
    {
        for(int j = 0; j<10; j++)
        {
            x[i][j] = (i * 10) + j + 1;
        }
    }



    ///////////////////////////////////////////////////
    //  print x  //////////////////////////////////////
    ///////////////////////////////////////////////////
    printf("/* print x[10][10] */\n");
    for(int i = 0; i< X_ROW ; i++)
    {
        for(int j = 0; j< X_ROW ; j++)
        {
            
            printf("x[%d][%d] : %d  ", i,j,x[i][j] );
            if(j==9) printf("\n");
        }
    }
    printf("\n");



    ///////////////////////////////////////////////////
    //  max pooling  //////////////////////////////////
    ///////////////////////////////////////////////////

    for(int i = 0; i< X_ROW/2; i++)
    {
        for(int j = 0; j<X_ROW/2; j++)
        {
            int max = x[2*i][2*j];
            
            for(int n = 0; n < POOLING_ROW; n++)
            {
                for(int m=0; m< POOLING_ROW; m++)
                {
                    if(max < x[2*i+n][2*j+m]){
                        max = x[2*i+n][2*j+m];
                    }
                }
            }

            x_pool[i][j] = max;
        }
        
    }



    ///////////////////////////////////////////////////
    //  print x_pool[i][j] ////////////////////////////
    ///////////////////////////////////////////////////
    printf("/* print x_pooling */\n");
    for(int i = 0; i < (X_ROW/2); i++)
    {
        for(int j = 0; j <(X_ROW/2); j++)
        {
            printf("x_pool[%d][%d] = %d  ", i,j,x_pool[i][j]);
            if(j==((X_ROW/2) -1)) printf("\n");
        }
        
    }


    return 0;
}