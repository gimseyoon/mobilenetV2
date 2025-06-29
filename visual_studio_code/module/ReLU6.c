#include <stdio.h>

#define ROW 14
#define COLUMN 14

void ReLU6(float *data);

int main()
{
    float data[384 * ROW * COLUMN] = {0};

    //////////////////////////////////////////////
    // Initialize input data /////////////////////
    //////////////////////////////////////////////
    for (int i = 0; i < 384 * ROW * COLUMN; i++) {
        data[i] = (float)i * 0.01f - 2.0f;  // 테스트용 : -2 ~ 4 정도 값 생성
    }

    ReLU6(data);

    //////////////////////////////////////////////
    // print result //////////////////////////////
    //////////////////////////////////////////////
    printf("\n=== ReLU6 result (first 1000 values) ===\n");
    for (int i = 0; i < 1000; i++) {
        printf("%f ", data[i]);
        if((i%10==0) && (i!=0)) printf("\n");
    }
    printf("\n");

    return 0;
}


void ReLU6(float *data)
{
    int CHANNEL = 384;

    for (int i = 0; i < CHANNEL; i++) {      /* channel */
        for (int j = 0; j < ROW; j++) {      /* row */
            for (int k = 0; k < COLUMN; k++) { /* column */
                int index = i * ROW * COLUMN + j * COLUMN + k;

                if (data[index] < 0)
                    data[index] = 0.0f;
                else if (data[index] > 6)
                    data[index] = 6.0f;
                // 0 <= data[index] <= 6 일 때는 그대로 유지
            }
        }
    }
}
