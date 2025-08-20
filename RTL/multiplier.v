// Module multiplier:
// receives 18-bit * [196] mul_in, multiplies each by a 17-bit mul_weight,
// rounds & shifts to 18-bit, and outputs as [196] mul_out.

`timescale 1ns / 1ps

module multiplier #(
    parameter IO_WIDTH      = 18,
    parameter ROW           = 14,
    parameter COLUMN        = 14,
    parameter PIXEL         = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH       = 17,
    parameter ADDR_CHANNEL  = $clog2(384),               // 8 (for CHANNEL = 384)
    parameter ADDR_WMEM     = $clog2(384 * 64),          // 15 (for 64*384 = 24576)
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
reg signed [IO_WIDTH*PIXEL-1:0] mul_in_q;
reg  signed [IO_WIDTH-1:0]  mul_out_reg [0:PIXEL-1];
wire signed [34:0]          mul_out_w   [0:PIXEL-1];

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
// local_state
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        mul_in_q <= 0;
    end else begin
        mul_in_q <= mul_in;
    end
end


///////////////////////////////////////////////////////
// Register outputs
///////////////////////////////////////////////////////
integer k;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        for (k = 0; k < PIXEL; k = k + 1)
            mul_out_reg[k] <= 0;
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
///////////////////////////////////////////////////////
genvar i;
generate
  for (i=0; i<PIXEL; i=i+1) begin : MULS
    // weight를 DSP 바로 앞에서 로컬화 (팬아웃=1)
    (* DONT_TOUCH = "true" *) reg signed [W_WIDTH-1:0] mul_weight_q;
    always @(posedge clk or negedge rst_n)
      if(!rst_n) mul_weight_q <= 0; else mul_weight_q <= mul_weight; 

    (* use_dsp = "yes" *)
    mult_gen_0 multiplier_0 (
      .CLK(clk),  // input wire CLK
      .A(mul_in_q[IO_WIDTH*(i+1)-1 : IO_WIDTH*i]),      // input wire [17 : 0] A
      .B(mul_weight_q),      // input wire [16 : 0] B
      .P(mul_out_w[i])      // output wire [34 : 0] P
    );

  end
endgenerate

endmodule