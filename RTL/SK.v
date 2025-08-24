`timescale 1ns / 1ps
module SK #(
    parameter IO_WIDTH = 18,
    parameter ROW      = 14,
    parameter COLUMN   = 14,
    parameter PIXEL    = ROW * COLUMN            // 196
)(
    input                                   clk,
    input                                   rst_n,
    input                                   skip_valid,                 // 1-cycle start pulse
    input       signed [IO_WIDTH*PIXEL-1:0] in_1,
    input       signed [IO_WIDTH*PIXEL-1:0] in_2,
    output      signed [IO_WIDTH*PIXEL-1:0] result,                     // 196개 패킹
    output reg                              result_save_valid,          // 1clk 펄스
    output reg                              skip_done                   // 64블록마다 1clk 펄스
);
    // -------------------------
    // 제어/카운터
    // -------------------------
    reg  busy;
    reg  skip_valid_d;
    wire start = skip_valid & ~skip_valid_d;

    // one-hot 쓰기 위치 (K): wr_sel는 현재 K, wr_sel_q는 K-1 (쓰기용)
    reg  [PIXEL-1:0] wr_sel;
    reg  [PIXEL-1:0] wr_sel_q;
    wire             last_write    = wr_sel[PIXEL-1]; // K==195?
    reg              last_write_q;

    // 블록(0..63) 카운터
    reg  [6:0] sk_cnt;
    reg        skip_done_pending;

    // -------------------------
    // 입력을 시프트로 순차 선택 (동적 part-select 제거)
    // -------------------------
    reg  [IO_WIDTH*PIXEL-1:0] in1_shift, in2_shift;
    wire signed [IO_WIDTH-1:0] in1_head = in1_shift[IO_WIDTH-1:0];
    wire signed [IO_WIDTH-1:0] in2_head = in2_shift[IO_WIDTH-1:0];

    // 합산 파이프(1클럭): 현재 합과 직전 합
    reg  signed [IO_WIDTH-1:0] sum_q;
    reg  signed [IO_WIDTH-1:0] sum_prev;

    // -------------------------
    // 결과 버퍼 & 패킹
    // -------------------------
    reg  signed [IO_WIDTH-1:0] result_reg [0:PIXEL-1];
    wire [IO_WIDTH*PIXEL-1:0]  packed_result;
    reg  [IO_WIDTH*PIXEL-1:0]  result_q;
    assign result = result_q;

    genvar g;
    generate
        for (g = 0; g < PIXEL; g = g + 1) begin : PACK_RESULT
            assign packed_result[IO_WIDTH*(g+1)-1 : IO_WIDTH*g] = result_reg[g];
        end
    endgenerate

    reg save_pending;
    integer k;

    // -------------------------
    // 메인 시퀀스
    // -------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            busy               <= 1'b0;
            skip_valid_d       <= 1'b0;

            wr_sel             <= {PIXEL{1'b0}};
            wr_sel_q           <= {PIXEL{1'b0}};
            last_write_q       <= 1'b0;

            in1_shift          <= {IO_WIDTH*PIXEL{1'b0}};
            in2_shift          <= {IO_WIDTH*PIXEL{1'b0}};
            sum_q              <= {IO_WIDTH{1'b0}};
            sum_prev           <= {IO_WIDTH{1'b0}};

            for (k = 0; k < PIXEL; k = k + 1) result_reg[k] <= {IO_WIDTH{1'b0}};

            result_q           <= {IO_WIDTH*PIXEL{1'b0}};
            result_save_valid  <= 1'b0;

            sk_cnt             <= 7'd0;
            skip_done          <= 1'b0;
            skip_done_pending  <= 1'b0;

            save_pending       <= 1'b0;
        end else begin
            skip_valid_d      <= skip_valid;
            result_save_valid <= 1'b0;
            skip_done         <= 1'b0;

            // ---------------------
            // 시작: 0번째 합을 "즉시" 계산 + 쉬프터는 한 칸 미리 당겨 놓음
            // ---------------------
            if (start && !busy) begin
                // 다음 사이클부터 head가 1번째 요소가 되도록 미리 1칸 쉬프트된 값 저장
                in1_shift <= { {IO_WIDTH{1'b0}}, in_1[IO_WIDTH*PIXEL-1:IO_WIDTH] };
                in2_shift <= { {IO_WIDTH{1'b0}}, in_2[IO_WIDTH*PIXEL-1:IO_WIDTH] };

                busy      <= 1'b1;

                // K=0부터 진행. 첫 사이클은 wr_sel_q=0 → 쓰기 없음
                wr_sel       <= {{(PIXEL-1){1'b0}}, 1'b1};
                wr_sel_q     <= {PIXEL{1'b0}};
                last_write_q <= 1'b0;

                // 0번째 요소의 합을 즉시 계산해 sum_q에 적재
                sum_q    <= $signed(in_1[IO_WIDTH-1:0]) + $signed(in_2[IO_WIDTH-1:0]);
                sum_prev <= {IO_WIDTH{1'b0}};
            end
            // ---------------------
            // 바쁠 때: K-1 쓰고, K 합 계산, 쉬프트/포인터 진행
            // ---------------------
            else if (busy) begin
                // 1) 직전 합을 K-1 위치에 기록 (정렬 보정)
                for (k = 0; k < PIXEL; k = k + 1)
                    if (wr_sel_q[k]) result_reg[k] <= sum_prev;

                // 2) 합 파이프: 현재 것을 prev로 넘기고 새 합 계산(현재 head 기준)
                sum_prev <= sum_q;
                sum_q    <= in1_head + in2_head;

                // 3) 입력 쉬프트 (다음 요소 준비)
                in1_shift <= { {IO_WIDTH{1'b0}}, in1_shift[IO_WIDTH*PIXEL-1:IO_WIDTH] };
                in2_shift <= { {IO_WIDTH{1'b0}}, in2_shift[IO_WIDTH*PIXEL-1:IO_WIDTH] };

                // 4) 쓰기 위치/마지막 플래그 1클럭 정렬
                wr_sel_q     <= wr_sel;
                last_write_q <= last_write;
                wr_sel       <= {wr_sel[PIXEL-2:0], 1'b0};

                // 5) 직전 사이클이 K=195 쓰기였다면 블록 종료
                if (last_write_q) begin
                    busy              <= 1'b0;
                    save_pending      <= 1'b1;              // 다음 사이클에 패킹 & valid
                    skip_done_pending <= (sk_cnt == 7'd63);

                    // 64블록 카운트
                    if (sk_cnt == 7'd63) sk_cnt <= 7'd0;
                    else                  sk_cnt <= sk_cnt + 1'b1;
                end
            end

            // 블록 완료: 패킹 및 유효 펄스
            if (save_pending) begin
                result_q          <= packed_result;   // 196개 패킹
                result_save_valid <= 1'b1;            // 1clk 펄스
                skip_done         <= skip_done_pending;
                skip_done_pending <= 1'b0;
                save_pending      <= 1'b0;
            end
        end
    end
endmodule
