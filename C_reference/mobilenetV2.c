#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//Parameter
#define EPS 1e-05   //epsilon for Batch Normalization
#define ROW 14      //Input Row
#define COLUMN 14   //Input Column
#define DK 3        //Depthwise Kernel Size
#define M 64        //Number of IN_Channel
#define N 384       //Number of Out_Channel
#define Z 1         //zero padding
#define S 1         //stride


//////////////////////////////////////////////////////////////////////////////////////////////////////////
// PATH of FILE  /////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

//INPUT FILE
#define LAYER_07_RESULT_FILE "C:/seyoon/mobilenetV2/value/result/provided_result/bin/mobilenetv2_result_layer07.bin"
#define INPUT_FILE "C:/seyoon/mobilenetV2/value/result/reference_result/layer_8/input.txt"
//OUTPUT FILE
#define OUTPUT_1_FILE "C:/seyoon/mobilenetV2/value/result/reference_result/layer_8/output_1.txt"
#define OUTPUT_1_BN_FILE "C:/seyoon/mobilenetV2/value/result/reference_result/layer_8/output_1_bn.txt"
#define OUTPUT_1_RELU_FILE "C:/seyoon/mobilenetV2/value/result/reference_result/layer_8/output_1_relu.txt"
#define OUTPUT_1_RELU_PADDED_FILE "C:/seyoon/mobilenetV2/value/result/reference_result/layer_8/output_1_padded.txt"
#define OUTPUT_2_FILE "C:/seyoon/mobilenetV2/value/result/reference_result/layer_8/output_2.txt"
#define OUTPUT_2_BN_FILE "C:/seyoon/mobilenetV2/value/result/reference_result/layer_8/output_2_bn.txt"
#define OUTPUT_2_RELU_FILE "C:/seyoon/mobilenetV2/value/result/reference_result/layer_8/output_2_relu.txt"
#define OUTPUT_3_FILE "C:/seyoon/mobilenetV2/value/result/reference_result/layer_8/output_3.txt"
#define OUTPUT_3_BN_FILE "C:/seyoon/mobilenetV2/value/result/reference_result/layer_8/output_3_bn.txt"
#define OUTPUT_3_SKIPCONNECTION_FILE "C:/seyoon/mobilenetV2/value/result/reference_result/layer_8/output_3_skip_connection.txt"
//WEIGHT FILE
#define PW_1_W_FILE "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_0_0_weight_24576.txt"
#define PW_2_W_FILE "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_2_weight_24576.txt"
#define DW_W_FILE "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_1_0_weight_3456.txt"
//PARAMETER FILE
#define MEAN_PATH_1 "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_0_1_running_mean_384.txt"
#define VAR_PATH_1 "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_0_1_running_var_384.txt"
#define BIAS_PATH_1 "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_0_1_bias_384.txt"
#define WEIGHT_PATH_1 "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_0_1_weight_384.txt"
#define MEAN_PATH_2 "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_1_1_running_mean_384.txt"
#define VAR_PATH_2 "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_1_1_running_var_384.txt"
#define BIAS_PATH_2 "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_1_1_bias_384.txt"
#define WEIGHT_PATH_2 "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_1_1_weight_384.txt"
#define MEAN_PATH_3 "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_3_running_mean_64.txt"
#define VAR_PATH_3 "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_3_running_var_64.txt"
#define BIAS_PATH_3 "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_3_bias_64.txt"
#define WEIGHT_PATH_3 "C:/seyoon/mobilenetV2/value/weight/layer_8/features_8_conv_3_weight_64.txt"


//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function Declare /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void CHECK_OUTPUT_WITH_TXT_FILE( const char *filename, const float *data, int channel);
void POINTWISE(int pw_num, float *in, float *out);
void DEPTHWISE(float *in, float *out);
void BATCH_NORMALIZATION(int bn_num, float *data);
void RELU6(float *data, int channel);
void SKIP_CONNECTION(float *input, float *output_3, float *output);


//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main Function /////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
    float input[M][ROW][COLUMN] = {0};      // input : result of layer_7
    float output_1[N][ROW][COLUMN]={0};     // output of 1st pointwise
    float output_2[N][ROW][COLUMN]={0};     // output of depthwise
    float output_3[M][ROW][COLUMN]={0};     // output of 2nd pointwise
    float output[M][ROW][COLUMN]={0};       // output : result of layer_8


    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // READ INPUT FILE ///////// /////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    FILE *fp_input = fopen(LAYER_07_RESULT_FILE, "rb");
    if(fp_input == NULL) {
        printf("Error: Cannot open input file.\n");
        return 1;
    }
    fread((void *)input, sizeof(float), M * ROW * COLUMN, fp_input); 
    fclose(fp_input);
    
    CHECK_OUTPUT_WITH_TXT_FILE(INPUT_FILE, &input[0][0][0],M);



    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1st pointwise convolution /////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    
    POINTWISE( 1, &input[0][0][0], &output_1[0][0][0] );
    CHECK_OUTPUT_WITH_TXT_FILE(OUTPUT_1_FILE, &output_1[0][0][0],N);

    BATCH_NORMALIZATION( 1, &output_1[0][0][0] );
    CHECK_OUTPUT_WITH_TXT_FILE(OUTPUT_1_BN_FILE, &output_1[0][0][0],N);
    
    RELU6( &output_1[0][0][0] , N);
    CHECK_OUTPUT_WITH_TXT_FILE(OUTPUT_1_RELU_FILE, &output_1[0][0][0],N);
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // depthwise convolution /////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////  
    
    DEPTHWISE( &output_1[0][0][0], &output_2[0][0][0]);
    CHECK_OUTPUT_WITH_TXT_FILE(OUTPUT_2_FILE, &output_2[0][0][0],N);

    BATCH_NORMALIZATION( 2, &output_2[0][0][0] );
    CHECK_OUTPUT_WITH_TXT_FILE(OUTPUT_2_BN_FILE, &output_2[0][0][0],N);

    RELU6( &output_2[0][0][0], N );
    CHECK_OUTPUT_WITH_TXT_FILE(OUTPUT_2_RELU_FILE, &output_2[0][0][0],N);



    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2nd pointwise convolution /////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    POINTWISE( 2, &output_2[0][0][0], &output_3[0][0][0] );
    CHECK_OUTPUT_WITH_TXT_FILE(OUTPUT_3_FILE, &output_3[0][0][0],M);

    BATCH_NORMALIZATION( 3, &output_3[0][0][0] );
    CHECK_OUTPUT_WITH_TXT_FILE(OUTPUT_3_BN_FILE, &output_3[0][0][0],M);
    
    SKIP_CONNECTION( &input[0][0][0], &output_3[0][0][0], &output[0][0][0] );
    CHECK_OUTPUT_WITH_TXT_FILE(OUTPUT_3_SKIPCONNECTION_FILE, &output[0][0][0],M);



    return 0;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTION //////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void CHECK_OUTPUT_WITH_TXT_FILE( const char *filename, const float *data, int channel)
{
    FILE *fp = fopen(filename, "wt");
    if (fp == NULL) {
        printf("error : CHECK_OUTPUT_WITH_TXT_FILE");
        return;
    }
    for (int c = 0; c < channel; c++) {
        for (int r = 0; r < ROW; r++) {
            for (int k = 0; k < COLUMN; k++) {
                /*fprintf(fp, "out[%d][%d][%d] : %f\n", c, r, k, data[c * ROW * COLUMN + r * COLUMN + k]);*/
                fprintf(fp, "%f\n", data[c * ROW * COLUMN + r * COLUMN + k]);
            }
        }
    }
    fclose(fp);
}
void POINTWISE(int pw_num, float *in, float *out)
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // SET (in_channel, out channel) /////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    int in_channel, out_channel;
    if(pw_num == 1) {
        in_channel = M;     //64
        out_channel = N;    //384
    }
    else if(pw_num == 2){
        in_channel = N;     //384   
        out_channel = M;    //64
    }
    float *weight = (float *)malloc(sizeof(float)*in_channel*out_channel);

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // READ weight_file //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    // 1st pointwise's weight
    if(pw_num == 1){        
        FILE *fp_weight = fopen(PW_1_W_FILE,"rt");
        if(fp_weight == NULL) {
            printf("Error: Cannot open weight_1 file.\n");
        }
        for (int i = 0; i < out_channel; i++) {
            for (int j = 0; j < in_channel; j++) {
                fscanf(fp_weight, "%f,", &weight[(in_channel*i) + j]);
            }
        }
        fclose(fp_weight);  
    }

    // 2nd pointwise's weight
    else if(pw_num == 2){   
        FILE *fp_weight = fopen(PW_2_W_FILE,"rt");
        if(fp_weight == NULL) {
            printf("Error: Cannot open weight_2 file.\n");
        }
        for (int i = 0; i < out_channel; i++) {
            for (int j = 0; j < in_channel; j++) {
                fscanf(fp_weight, "%f,", &weight[(in_channel*i) + j]);
            }
        }
        fclose(fp_weight);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // pointwise convolution /////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    for(int i = 0; i < out_channel; i++) {
        for(int k = 0; k < ROW; k++) {
            for(int p = 0; p < COLUMN; p++) {
                for(int j = 0; j < in_channel; j++) {
                    out[(i * ROW * COLUMN) + (k * COLUMN) + p] += 
                    in[(j * ROW * COLUMN) + (k * COLUMN) + p] * 
                    weight[(i * in_channel) + j];
                }
            }
        }
    }

    free(weight);
}
void DEPTHWISE(float *in, float *out)
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //declare/////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////// 

    float output_1_padded[N][ROW + 2*Z][COLUMN + 2*Z] = {0};
    float weight[N][DK][DK] = {0};



    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // READ WEIGHT FILE //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    FILE *fp_weight = fopen(DW_W_FILE,"rt");
    if(fp_weight == NULL) {
        printf("Error: Cannot open weight_2 file.\n");
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < DK; j++) {
            for(int k = 0; k < DK; k++){
                fscanf(fp_weight, "%f,", &weight[i][j][k]);
            }
        }
    }
    fclose(fp_weight);  

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // padding part //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////     

    for(int i=0; i < N; i++){
        for(int j = 0; j < ROW; j++){
            for(int k = 0; k < COLUMN; k++){
                output_1_padded[i][ j + Z ][ k + Z ] = in[ (i*ROW*COLUMN) + (j*COLUMN) + (k) ];
            }
        }
    }
    
    FILE *fp_pad = fopen(OUTPUT_1_RELU_PADDED_FILE, "wt");
    if (fp_pad == NULL) {
        printf("Error: Cannot open padded input output file.\n");
        return;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < ROW + 2*Z; j++) {
            for (int k = 0; k < COLUMN + 2*Z; k++) {
                fprintf(fp_pad, "%.6f\n", output_1_padded[i][j][k]);
            }
        }
    }
    fclose(fp_pad);

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // compute part //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////     

    for(int i=0; i < N; i++){     // Number of multiplication operations : N * DF * DF * DK * DK
        for(int j = 0; j < ROW; j++) {
            for(int k = 0; k < COLUMN; k++) {
                for(int p = 0; p < DK; p++) {
                    for(int u = 0; u < DK; u++){
                        out[ (i*ROW*COLUMN) + (j*COLUMN) + (k)] += output_1_padded[i][ j +p ][ k+u ] * weight[i][p][u];
                    }
                }     
            }  
        }
    }
}
void BATCH_NORMALIZATION(int bn_num, float *data)
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // declare ///////////////////////////////////////////////////////////////////////////////////////////   
    //////////////////////////////////////////////////////////////////////////////////////////////////////   
    const char *mean_path;
    const char *var_path;
    const char *bias_path;
    const char *weight_path;
    int CHANNEL;



    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // SET 'CHANNEL' and 'file path'  /////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////    
    if (bn_num == 1) {
        CHANNEL = N;
        mean_path   = MEAN_PATH_1;
        var_path    = VAR_PATH_1;
        bias_path   = BIAS_PATH_1;
        weight_path = WEIGHT_PATH_1;
    } 
    else if (bn_num == 2) {
        CHANNEL = N;
        mean_path   = MEAN_PATH_2;
        var_path    = VAR_PATH_2;
        bias_path   = BIAS_PATH_2;
        weight_path = WEIGHT_PATH_2;
    } 
    else if (bn_num == 3) {
        CHANNEL = M;
        mean_path   = MEAN_PATH_3;
        var_path    = VAR_PATH_3;
        bias_path   = BIAS_PATH_3;
        weight_path = WEIGHT_PATH_3;
    }



    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // malloc for parameter //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    float *mean = (float*)malloc(CHANNEL * sizeof(float));
    float *var = (float*)malloc(CHANNEL * sizeof(float));
    float *bias = (float*)malloc(CHANNEL * sizeof(float));
    float *weight = (float*)malloc(CHANNEL * sizeof(float));



    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // READ parameter file ///////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    FILE *fp_mean = fopen(mean_path, "rt");
    FILE *fp_var = fopen(var_path, "rt");
    FILE *fp_bias = fopen(bias_path, "rt");
    FILE *fp_weight = fopen(weight_path, "rt");

    if (!fp_mean || !fp_var || !fp_bias || !fp_weight) {
        printf("Error opening parameter files.\n");
        exit(1);
    }

    for (int i = 0; i < CHANNEL; i++) {
        fscanf(fp_mean, "%f,", &mean[i]);
        fscanf(fp_var, "%f,", &var[i]);
        fscanf(fp_bias, "%f,", &bias[i]);
        fscanf(fp_weight, "%f,", &weight[i]);
    }

    fclose(fp_mean);
    fclose(fp_var);
    fclose(fp_bias);
    fclose(fp_weight);



    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // Batch Normalization ///////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int i = 0; i < CHANNEL; i++) {         /* channel */
        for (int j = 0; j < ROW; j++) {         /* row */
            for (int k = 0; k < COLUMN; k++) {  /* column */
                int index = (i * ROW * COLUMN) + (j * COLUMN) + (k);
                data[index] = ((data[index] - mean[i]) / sqrt(var[i] + EPS)) * weight[i] + bias[i];
            }
        }
    }



    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // free(parameter) ///////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    free(mean);
    free(var);
    free(bias);
    free(weight);
}
void RELU6(float *data, int channel) {
    for (int i = 0; i < channel; i++) {
        for (int j = 0; j < ROW; j++) {
            for (int k = 0; k < COLUMN; k++) {
                int index = (i * ROW * COLUMN) + (j * COLUMN) + k;
                if (data[index] < 0) data[index] = 0.0f;
                else if (data[index] > 6) data[index] = 6.0f;
            }
        }
    }
}
void SKIP_CONNECTION(float *data_1, float *data_2, float *out)
{
    for (int i = 0; i < M; i++) {
        for(int j = 0; j < ROW; j++){
            for(int k = 0; k < COLUMN; k++){
                out[ (i*ROW*COLUMN) + (j*COLUMN) + (k)] = data_1[(i*ROW*COLUMN) + (j*COLUMN) + (k)] 
                                                        + data_2[(i*ROW*COLUMN) + (j*COLUMN) + (k)];
            }
        }
    }
}
