#include <stdio.h>
#define INPUT_SIZE 5
#define FILTER_SIZE 3

int main()
{
    int x[5] = {7,2,3,3,8};
    int w[3] = {1,0,-1};
    int y[3] = {0, };

    for(int i = 0; i < INPUT_SIZE; i++){
        for(int j = 0; j < FILTER_SIZE;j++){
            y[i] += x[i+j] * w[j] ; 
        }
    }

    for(int k = 0; k < (INPUT_SIZE - FILTER_SIZE + 1) ; k++){
        printf("y[%d] = %d\n", k, y[k] );
    
    }

    return 0;
} 
