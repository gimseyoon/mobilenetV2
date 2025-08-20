// Module multiplier:
// receives 18-bit * [196] mul_in, multiplies each by a 17-bit mul_weight,
// rounds & shifts to 18-bit, and outputs as [196] mul_out.

`timescale 1ns / 1ps

module multiplier #(
    parameter IO_WIDTH       = 18,
    parameter ROW            = 14,
    parameter COLUMN         = 14,
    parameter PIXEL          = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH        = 17,
    parameter ADDR_CHANNEL   = $clog2(384),               // 8 (for CHANNEL = 384)
    parameter ADDR_WMEM      = $clog2(384 * 64),          // 15 (for 64*384 = 24576)
    parameter integer R_SHIFT = 16
)(
    input                                       clk,
    input                                       rst_n,
    input                                       pw_1_done,
    input                                       dw_done,
    input                                       pw_2_done,
    input  signed [IO_WIDTH*PIXEL-1:0]          mul_in,      // [3528-1:0]
    input  signed [W_WIDTH-1:0]                 mul_weight,  // [16:0]
    output signed [IO_WIDTH*PIXEL-1:0]          mul_out      // [3528-1:0]
);

///////////////////////////////////////////////////////
// Internal regs/wires
///////////////////////////////////////////////////////
reg  signed [IO_WIDTH*PIXEL-1:0] mul_in_q;

// 행 단위 가중치 복제 (fan-out 196 -> 14)
reg  signed [W_WIDTH-1:0]  mul_weight_row_q [0:ROW-1];

reg  signed [IO_WIDTH-1:0] mul_out_reg [0:PIXEL-1];
wire signed [34:0]         mul_out_w   [0:PIXEL-1];

///////////////////////////////////////////////////////
// Rounding + right-shift to signed 18-bit
///////////////////////////////////////////////////////
function signed [17:0] round_shift_signed18;
    input signed [34:0] x;
    reg        sign;
    reg [34:0] mag;
    reg [35:0] mag_rnd;
    reg [17:0] q_u;
begin
    sign    = x[34];
    mag     = sign ? (~x + 35'd1) : x;              // |x|
    mag_rnd = {1'b0, mag} + (36'd1 << (R_SHIFT-1)); // +0.5LSB
    q_u     = mag_rnd[35:R_SHIFT];                  // >> R
    round_shift_signed18 = sign ? -$signed(q_u) : $signed(q_u);
end
endfunction

///////////////////////////////////////////////////////
// Pipeline: mul_in (1clk)
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) mul_in_q <= {IO_WIDTH*PIXEL{1'b0}};
    else        mul_in_q <= mul_in;
end

///////////////////////////////////////////////////////
// Pipeline: mul_weight (행 단위로 1clk 복제)
///////////////////////////////////////////////////////
integer r;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        for (r=0; r<ROW; r=r+1) mul_weight_row_q[r] <= {W_WIDTH{1'b0}};
    end else begin
        // 한 사이클에 14개 행 레지스터 모두 동일한 mul_weight 캡처
        for (r=0; r<ROW; r=r+1) mul_weight_row_q[r] <= mul_weight;
    end
end
// 참고: 필요시 아래와 같이 도와줄 수도 있음 (자동 복제 유도)
// (* max_fanout = 16 *) wire signed [W_WIDTH-1:0] mul_weight_buf = mul_weight;

///////////////////////////////////////////////////////
// Register outputs
///////////////////////////////////////////////////////
integer k;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        for (k = 0; k < PIXEL; k = k + 1)
            mul_out_reg[k] <= {IO_WIDTH{1'b0}};
    end else begin
        for (k = 0; k < PIXEL; k = k + 1)
            mul_out_reg[k] <= round_shift_signed18(mul_out_w[k]);
    end
end

///////////////////////////////////////////////////////
// Pack register array -> bus
///////////////////////////////////////////////////////
genvar m;
generate
    for (m = 0; m < PIXEL; m = m + 1) begin : OUTPUT_PACK
        assign mul_out[IO_WIDTH*(m+1)-1 : IO_WIDTH*m] = mul_out_reg[m];
    end
endgenerate

///////////////////////////////////////////////////////
// MULTIPLIER instantiation (per pixel)
//   - B 입력은 "행 레지스터"에서만 가져옴 (fan-out ~14)
//   - A 입력은 mul_in_q에서 픽셀별 슬라이스
///////////////////////////////////////////////////////
genvar ry, cx;
generate
    for (ry = 0; ry < ROW;    ry = ry + 1) begin : ROWS
    for (cx = 0; cx < COLUMN; cx = cx + 1) begin : COLS
        localparam integer IDX = ry*COLUMN + cx;

        (* use_dsp = "yes" *)
        mult_gen_0 u_mul (
            .CLK (clk),
            .A   ($signed(mul_in_q[IO_WIDTH*(IDX+1)-1 : IO_WIDTH*IDX])),
            .B   (mul_weight_row_q[ry]),   // 행 단위 1clk 레지스터
            .P   (mul_out_w[IDX])
        );
    end
    end
endgenerate

endmodule
