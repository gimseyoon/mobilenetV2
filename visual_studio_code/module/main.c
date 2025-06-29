#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 1e-05
#define ROW 14
#define COLUMN 14
#define M 64
#define N 384

void POINTWISE(int pw_num, float *input, float *output_1);
void DEPTHWISE(float *output_1, float *output_2);
void BATCH_NORMALIZATION(int bn_num, float *data);
void RELU6(float *data);
void SKIP_CONNECTION(float *input_data, float *output_3);


int main()
{
    float input[M][ROW][COLUMN] = {0};
    float output_1[N][ROW][COLUMN]={0};    //output of 1st pointwise
    float output_2[N][ROW][COLUMN]={0};    //output of depthwise
    float output_3[M][ROW][COLUMN]={0};     //output of 2nd pointwise

    //////////////////////////////////
    // READ input file ///////////////
    //////////////////////////////////
    FILE *fp_input = fopen("C:/seyoon/parameter/data/bin","rb");
    if(fp_input == NULL) {
        printf("Error: Cannot open input file.\n");
        return 1;
    }
    fread((void *)input, sizeof(float), M * ROW * COLUMN, fp_input); 
    fclose(fp_input);

    //////////////////////////////////
    // 1st pointwise convolution /////
    //////////////////////////////////
    POINTWISE( 1, input, output_1 );
    BATCH_NORMALIZATION( 1, output_1 );
    RELU6( output_1 );

    //////////////////////////////////
    // depthwise convolution /////////
    //////////////////////////////////  
    DEPTHWISE( output_1, output_2);
    BATCH_NORMALIZATION( 2, output_2 );
    RELU6( output_2 );
    
    //////////////////////////////////
    // 2nd pointwise convolution /////
    //////////////////////////////////
    POINTWISE( 2, output_2, output_3 );
    BATCH_NORMALIZATION( 1, output_3 );
    SKIP_CONNECTION( input, output_3 );


    return 0;
}