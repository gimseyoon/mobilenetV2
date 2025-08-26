`timescale 1ns / 1ps

module multiplier #(
    parameter IO_WIDTH        = 18,
    parameter ROW             = 14,
    parameter COLUMN          = 14,
    parameter PIXEL           = ROW * COLUMN,   // 196
    parameter W_WIDTH         = 17,
    parameter integer R_SHIFT = 16,
    parameter integer SUB     = 2
)(
    input                               clk,
    input                               rst_n,
    input  signed [IO_WIDTH*PIXEL-1:0]  mul_in,
    input  signed [W_WIDTH-1:0]         mul_weight,
    output signed [IO_WIDTH*PIXEL-1:0]  mul_out
);

///////////////////////////////////////////////////////
reg signed [IO_WIDTH*PIXEL-1:0] mul_in_q;
always @(posedge clk or negedge rst_n)
    if (!rst_n) mul_in_q <= {IO_WIDTH*PIXEL{1'b0}};
    else        mul_in_q <= mul_in;

///////////////////////////////////////////////////////
// weight 
///////////////////////////////////////////////////////
(* DONT_TOUCH = "true" *)
reg signed [W_WIDTH-1:0] w_dup [0:ROW*SUB-1];

integer rr, gg;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        for (rr=0; rr<ROW*SUB; rr=rr+1) w_dup[rr] <= {W_WIDTH{1'b0}};
    end else begin
        // ???? weight?? ?????? ????? ????(?????????)
        for (rr=0; rr<ROW; rr=rr+1)
            for (gg=0; gg<SUB; gg=gg+1)
                w_dup[rr*SUB + gg] <= mul_weight;
    end
end


///////////////////////////////////////////////////////
function signed [17:0] round_shift_signed18;
    input signed [34:0] x;
    reg sign; reg [34:0] mag; reg [35:0] mag_rnd; reg [17:0] q_u;
begin
    sign    = x[34];
    mag     = sign ? (~x + 35'd1) : x;
    mag_rnd = {1'b0, mag} + (36'd1 << (R_SHIFT-1));
    q_u     = mag_rnd[35:R_SHIFT];
    round_shift_signed18 = sign ? -$signed(q_u) : $signed(q_u);
end
endfunction

reg  signed [IO_WIDTH-1:0] mul_out_reg [0:PIXEL-1];
wire signed [34:0]         mul_out_w   [0:PIXEL-1];

integer k;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        for (k=0; k<PIXEL; k=k+1) mul_out_reg[k] <= {IO_WIDTH{1'b0}};
    end else begin
        for (k=0; k<PIXEL; k=k+1) mul_out_reg[k] <= round_shift_signed18(mul_out_w[k]);
    end
end

genvar m;
generate
    for (m=0; m<PIXEL; m=m+1) begin : PACK
        assign mul_out[IO_WIDTH*(m+1)-1 : IO_WIDTH*m] = mul_out_reg[m];
    end
endgenerate


///////////////////////////////////////////////////////
genvar ry, cx;
generate
    for (ry=0; ry<ROW; ry=ry+1) begin : ROWS
        for (cx=0; cx<COLUMN; cx=cx+1) begin : COLS
            localparam integer IDX  = ry*COLUMN + cx;
          
            localparam integer GIDX = (cx*SUB)/COLUMN; // 0..SUB-1

            (* use_dsp = "yes" *)
            mult_gen_0 u_mul (
                .CLK (clk),
                .A   ($signed(mul_in_q[IO_WIDTH*(IDX+1)-1 : IO_WIDTH*IDX])),
                .B   (w_dup[ry*SUB + GIDX]),
                .P   (mul_out_w[IDX])
            );
        end
    end
endgenerate

endmodule