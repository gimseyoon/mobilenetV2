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
    output      signed [IO_WIDTH*PIXEL-1:0] result,                     // 완성된 프레임만 출력
    output reg                              result_save_valid,          // 1클럭 펄스
    output reg                              skip_done                   // 64프레임 완료 펄스
);
    // -------------------------
    // 상태/포인터/카운터
    // -------------------------
    reg        busy;
    reg  [7:0] idx;                        // 0..195
    reg        skip_valid_d;
    wire       start = skip_valid & ~skip_valid_d;

    reg  [PIXEL-1:0] wr_sel;               // one-hot write 포인터
    wire             last_write = wr_sel[PIXEL-1];

    reg  [6:0] sk_cnt;                     // 프레임 수(0..63)
    reg        skip_done_pending;          // result_save_valid와 동기화용

    // -------------------------
    // 입력 래치 & 선택
    // -------------------------
    reg  signed [IO_WIDTH*PIXEL-1:0] in_1_q, in_2_q;
    wire signed [IO_WIDTH-1:0] in1_sel = in_1_q[IO_WIDTH*idx +: IO_WIDTH];
    wire signed [IO_WIDTH-1:0] in2_sel = in_2_q[IO_WIDTH*idx +: IO_WIDTH];
    wire signed [IO_WIDTH-1:0] sum_w   = in1_sel + in2_sel;   // 1클럭 당김(즉시 기록)

    // -------------------------
    // 결과 레지스터 뱅크(픽셀별 1회 쓰기 후 유지)
    // -------------------------
    reg  signed [IO_WIDTH-1:0] result_reg [0:PIXEL-1];

    // -------------------------
    // 결과 버스(완성 시점에 한 번에 래치)
    // -------------------------
    wire [IO_WIDTH*PIXEL-1:0] packed_result;
    reg  [IO_WIDTH*PIXEL-1:0] result_q;
    assign result = result_q;

    genvar g;
    generate
        for (g = 0; g < PIXEL; g = g + 1) begin : PACK_RESULT
            assign packed_result[IO_WIDTH*(g+1)-1 : IO_WIDTH*g] = result_reg[g];
        end
    endgenerate

    // 마지막 쓰기 직후 다음 클럭에 버스를 내보내기 위한 플래그
    reg save_pending;

    integer k;

    // -------------------------
    // 메인 시퀀셜
    // -------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            busy               <= 1'b0;
            idx                <= 8'd0;
            wr_sel             <= {PIXEL{1'b0}};
            skip_valid_d       <= 1'b0;
            result_save_valid  <= 1'b0;
            skip_done          <= 1'b0;
            sk_cnt             <= 7'd0;
            skip_done_pending  <= 1'b0;
            in_1_q             <= {IO_WIDTH*PIXEL{1'b0}};
            in_2_q             <= {IO_WIDTH*PIXEL{1'b0}};
            save_pending       <= 1'b0;
            result_q           <= {IO_WIDTH*PIXEL{1'b0}};
            for (k = 0; k < PIXEL; k = k + 1) result_reg[k] <= {IO_WIDTH{1'b0}};
        end else begin
            skip_valid_d      <= skip_valid;
            result_save_valid <= 1'b0;
            skip_done         <= 1'b0;

            // 프레임 시작: 입력 버스 래치 & 포인터 초기화
            if (start && !busy) begin
                in_1_q <= in_1;
                in_2_q <= in_2;
                busy   <= 1'b1;
                idx    <= 8'd0;
                wr_sel <= {{(PIXEL-1){1'b0}}, 1'b1};  // LSB=1
            end
            // 프레임 진행
            else if (busy) begin
                // 현재 포인터가 가리키는 위치에 **즉시 기록**(1클럭 당김)
                for (k = 0; k < PIXEL; k = k + 1)
                    if (wr_sel[k]) result_reg[k] <= sum_w;

                // 마지막 픽셀을 지금 썼으면 다음 클럭에 결과버스 내보내도록 예약
                if (last_write) begin
                    busy              <= 1'b0;
                    wr_sel            <= {PIXEL{1'b0}};
                    save_pending      <= 1'b1;              // 다음 클럭에 pack & valid
                    skip_done_pending <= (sk_cnt == 7'd63); // 64번째였는지 저장

                    // 64프레임 카운트
                    if (sk_cnt == 7'd63) sk_cnt <= 7'd0;
                    else                  sk_cnt <= sk_cnt + 1'b1;
                end else begin
                    // 다음 픽셀로 진행 (진짜 one-hot: 0을 쉬프트 인)
                    wr_sel <= {wr_sel[PIXEL-2:0], 1'b0};
                    idx    <= idx + 1'b1;
                end
            end

            // 프레임 완성 → 다음 클럭에 한 줄로 묶어 출력 & valid 펄스
            if (save_pending) begin
                result_q          <= packed_result;     // 196개 한 번에 출력
                result_save_valid <= 1'b1;              // 1클럭 펄스
                skip_done         <= skip_done_pending; // 64번째 완료면 펄스
                skip_done_pending <= 1'b0;
                save_pending      <= 1'b0;
            end
        end
    end
endmodule
