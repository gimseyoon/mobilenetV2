#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 1e-05
#define ROW 14
#define COLUMN 14

void BN(int bn_num, float *data);

int main()
{
    float data1[384 * ROW * COLUMN] = {0};
    float data2[384 * ROW * COLUMN] = {0};
    float data3[64 * ROW * COLUMN] = {0};

    // data1 초기화 : 모든 값 1.0
    for (int i = 0; i < 384 * ROW * COLUMN; i++) {
        data1[i] = 1.0f;
    }

    // data2 초기화 : 인덱스 기반 증가값
    for (int i = 0; i < 384 * ROW * COLUMN; i++) {
        data2[i] = (float)i * 0.01f;
    }

    // data3 초기화 : 모두 0.5
    for (int i = 0; i < 64 * ROW * COLUMN; i++) {
        data3[i] = 0.5f;
    }

    BN(1, data1);
    BN(2, data2);
    BN(3, data3);

    //////////////////////////////////////////////
    // print result //////////////////////////////
    //////////////////////////////////////////////
    printf("\n=== data1 result (first 10 values) ===\n");
    for (int i = 0; i < 10; i++) {
        printf("%f ", data1[i]);
    }
    printf("\n");

    printf("\n=== data2 result (first 10 values) ===\n");
    for (int i = 0; i < 10; i++) {
        printf("%f ", data2[i]);
    }
    printf("\n");

    printf("\n=== data3 result (first 10 values) ===\n");
    for (int i = 0; i < 10; i++) {
        printf("%f ", data3[i]);
    }
    printf("\n");


    return 0;
}


void BN(int bn_num, float *data)
{
    ///////////////////////////////////////////////////
    // declare ////////////////////////////////////////
    ///////////////////////////////////////////////////     
    const char *mean_path;
    const char *var_path;
    const char *bias_path;
    const char *weight_path;
    int CHANNEL;



    ///////////////////////////////////////////////////
    // SET 'CHANNEL NUMBER' and 'file path'  //////////
    ///////////////////////////////////////////////////    
    if (bn_num == 1) {
        CHANNEL = 384;
        mean_path   = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/weight/features_8_conv_0_1_running_mean_384.txt";
        var_path    = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/weight/features_8_conv_0_1_running_var_384.txt";
        bias_path   = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/weight/features_8_conv_0_1_bias_384.txt";
        weight_path = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/weight/features_8_conv_0_1_weight_384.txt";
    } 
    else if (bn_num == 2) {
        CHANNEL = 384;
        mean_path   = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/weight/features_8_conv_1_1_running_mean_384.txt";
        var_path    = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/weight/features_8_conv_1_1_running_var_384.txt";
        bias_path   = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/weight/features_8_conv_1_1_bias_384.txt";
        weight_path = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/weight/features_8_conv_1_1_weight_384.txt";
    } 
    else if (bn_num == 3) {
        CHANNEL = 64;
        mean_path   = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/weight/features_8_conv_3_running_mean_64.txt";
        var_path    = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/weight/features_8_conv_3_running_var_64.txt";
        bias_path   = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/weight/features_8_conv_3_bias_64.txt";
        weight_path = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/weight/features_8_conv_3_weight_64.txt";
    }
    else {
        printf("Invalid bn_num: %d\n", bn_num);
        return;
    }


    ///////////////////////////////////////////////////
    // malloc for parameter ///////////////////////////
    ///////////////////////////////////////////////////
    float *mean = (float*)malloc(CHANNEL * sizeof(float));
    float *var = (float*)malloc(CHANNEL * sizeof(float));
    float *bias = (float*)malloc(CHANNEL * sizeof(float));
    float *weight = (float*)malloc(CHANNEL * sizeof(float));



    ///////////////////////////////////////////////////
    // READ parameter file ////////////////////////////
    ///////////////////////////////////////////////////
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



    ///////////////////////////////////////////////////
    // Batch Normalization ////////////////////////////
    ///////////////////////////////////////////////////
    for (int i = 0; i < CHANNEL; i++) {         /* channel */
        for (int j = 0; j < ROW; j++) {         /* row */
            for (int k = 0; k < COLUMN; k++) {  /* column */
                int index = i * ROW * COLUMN + j * COLUMN + k;
                data[index] = ((data[index] - mean[i]) / sqrt(var[i] + EPS)) * weight[i] + bias[i];
            }
        }
    }



    ///////////////////////////////////////////////////
    // free(parameter) ////////////////////////////////
    ///////////////////////////////////////////////////
    free(mean);
    free(var);
    free(bias);
    free(weight);
}
