 /*pointwise convolution 
 - input: [10][3][3]
 - output: [5][3][3]
 - weight: [5][10]

 + 과제: layer 08 <- layer 07의 결과를 입력으로
   float 형 부동소수점 결과라 좀 다를 수 있음
*/

#include <stdio.h>

#define DF 3
#define M 10
#define N 5



int main()
{
    int in[M][DF][DF] = {0};
    int w[N][M] = {0};
    int out[N][DF][DF] = {0};

    for(int i=0; i < M ; i++)
    {
        for(int j=0; j < DF; j++)
        {
            for(int k=0; k < DF; k++)
            {
                in[i][j][k] = 1;
            }
        }
    }

    for(int i=0; i < N ; i++)
    {
        for(int j=0; j < M; j++)
        {
            w[i][j]=1;
        }
    }


    for(int i=0; i < N ; i++){      // Number of multiplication operations : N * M * DF * DF
        for(int j=0; j < M; j++){
            for(int k=0; k < DF; k++){
                for(int p = 0; p < DF; p++){
                    out[i][k][p] += in[j][k][p] * w[i][j];
                }
            }
        }
    }

    for(int i=0; i < N ; i++){
        for(int j=0; j < DF; j++){
            for(int k=0; k < DF; k++){
                printf("y[%d][%d][%d] = %d\n",i,j,k, out[i][j][k]);
            }
        }
    }


    return 0;
}