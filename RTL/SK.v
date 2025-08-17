`timescale 1ns / 1ps

module SK #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN  // 196
)(
    input                                   clk,
    input                                   rst_n,
    input                                   skip_valid,                 // 1-cycle start pulse
    input       signed [IO_WIDTH*PIXEL-1:0] in_1,
    input       signed [IO_WIDTH*PIXEL-1:0] in_2,
    output      signed [IO_WIDTH*PIXEL-1:0] result,
    output reg                              result_save_valid,
    output reg                              skip_done
);


///////////////////////////////////////////////////////

reg [6:0] sk_cnt;


reg       busy;
reg [7:0] idx;         // 0..195
reg       skip_valid_d;
wire      start;

reg  signed [IO_WIDTH*PIXEL-1:0] in_1_q;
reg  signed [IO_WIDTH*PIXEL-1:0] in_2_q;
wire signed [IO_WIDTH-1:0]       in_1_array;
wire signed [IO_WIDTH-1:0]       in_2_array;
//(* use_dsp = "yes" *) 
wire signed [IO_WIDTH-1:0]       sum;
reg  signed [IO_WIDTH-1:0]       result_reg [0:PIXEL-1];


///////////////////////////////////////////////////////
assign start = skip_valid & ~skip_valid_d;
assign in_1_array = in_1_q[IO_WIDTH*idx +: IO_WIDTH];
assign in_2_array = in_2_q[IO_WIDTH*idx +: IO_WIDTH];
assign sum  = in_1_array + in_2_array;
/////////////////////////////////////////////////////// 
// Control: start, busy, index
///////////////////////////////////////////////////////
integer i;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        sk_cnt             <= 0;
        skip_valid_d       <= 0;
        busy               <= 0;
        idx                <= 0;
        result_save_valid  <= 0;
        in_1_q             <= 0;
        in_2_q             <= 0;
        skip_done          <= 0;
        for (i = 0; i < PIXEL; i = i + 1) begin
            result_reg[i] <= 0;
        end
    end 
    else begin
        skip_valid_d      <= skip_valid;
        result_save_valid <= 1'b0; // 기본 0, 완료 시 1로 1클럭

        // 시작: 입력 버스 래치 후 처리 시작
        if (start && !busy) begin
            in_1_q <= in_1;
            in_2_q <= in_2;
            busy   <= 1'b1;
            idx    <= 0;
            for (i = 0; i < PIXEL; i = i + 1) begin
                result_reg[i] <= 0;
            end
        end
        // 진행 중: 한 픽셀씩 처리
        else if (busy) begin
            // 현재 idx 위치에 결과 기록 (write-first)
            result_reg[idx] <= sum[IO_WIDTH-1:0];

            if (idx == 195) begin
                busy              <= 1'b0;      // 완료
                result_save_valid <= 1'b1;      // 결과 완성 신호 1클럭
                sk_cnt            <= sk_cnt + 1;
            end 
            else begin
                idx <= idx + 1'b1;
            end
        end
        else begin
            idx <= 0;
        end
        
        if(sk_cnt == 64 && idx == 195) begin
            skip_done <= 1;
        end
        else begin
            skip_done <= 0;
        end
    end
end


/////////////////////////////////////////////////////// 
// Pack result_reg[] -> result bus
///////////////////////////////////////////////////////
genvar m;
generate
    for (m = 0; m < PIXEL; m = m + 1) begin : PACK_RESULT
        assign result[IO_WIDTH*(m+1)-1 : IO_WIDTH*m] = result_reg[m];
    end
endgenerate

endmodule
