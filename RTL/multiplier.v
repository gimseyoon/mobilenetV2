// Module multiplier: 
// receives 20 bit * [196] mul_in, multiplies each by a 17 bit mul_weight,
// and outputs the result as 20bit [196] mul_out.

`timescale 1ns / 1ps

module multiplier #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),         // 8 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64),       // 15 (for 64*384 = 24576)
    parameter integer R_SHIFT = 16
)(
    input                                       clk,
    input                                       rst_n,
    input                                       pw_1_done,
    input                                       dw_done,
    input                                       pw_2_done,
    input   signed [IO_WIDTH * PIXEL - 1 : 0]   mul_in,  // [3528-1 : 0], 3528 bit
    input   signed [W_WIDTH - 1 : 0]            mul_weight,   // [17-1:0], 17 bit
    output  signed [IO_WIDTH * PIXEL - 1 : 0]   mul_out  // [3528-1 : 0], 3528 bit
    );
/////////////////////////////////////////////////////////////////////////

    reg signed [IO_WIDTH-1 :0] mul_out_reg [0 : PIXEL-1]; // [384][14][14]
    wire signed [34:0] mul_out_w [0 : PIXEL-1]; //[196-1 :0]
    
    
/////////////////////////////////////////////////////////////////////////
        
    function signed [17:0] round_shift_signed18;
        input signed [34:0] x;
        reg        sign;
        reg [34:0] mag;
        reg [35:0] mag_rnd;
        reg [17:0] q_u;
    begin
        sign    = x[34];
        mag     = sign ? (~x + 35'd1) : x;                  // |x|
        mag_rnd = {1'b0, mag} + (36'd1 << (R_SHIFT-1));     // +0.5LSB
        q_u     = mag_rnd[35:R_SHIFT];                      // >> R
        round_shift_signed18 = sign ? -$signed(q_u) : $signed(q_u);
    end
    endfunction
    
/////////////////////////////////////////////////////////////////////////
// concat : mul_out_reg <= { mul_out_w[k][36] , mul_out_w[k][34 :16] }

    integer k;
    always@(posedge clk or negedge rst_n) begin
        if(!rst_n) begin
            for(k=0; k < PIXEL; k = k+1) begin
                mul_out_reg[k] <= 0;
            end
        end
        else begin
            if(pw_1_done || dw_done || pw_2_done ) begin
                for(k=0; k < PIXEL; k = k+1) begin
                    mul_out_reg[k] <= 0;
                end
            end
            else begin
                for(k=0; k < PIXEL; k = k+1) begin
                    mul_out_reg[k] <= round_shift_signed18(mul_out_w[k]);
                end
            end

        end //else
    end //always
    
    
////////////////////////////////////////////////////////////////////////
// assign [3528-1 :0] mul_out = [18-1:0] mul_out_reg [196-1 :0]

    genvar m;
    generate
        for (m = 0; m < PIXEL; m = m + 1) begin : OUTPUT_PACK
            assign mul_out[IO_WIDTH*(m+1)-1 : IO_WIDTH*m] = mul_out_reg[m];
        end
    endgenerate
    
    
    
    
    
//////////////////////////////////////////////////////////////////////
// MULTIPLIER instantiation

    genvar i;
    generate 
      for(i = 0; i < PIXEL; i = i + 1) begin
        (* use_dsp = "yes" *) mult_gen_0 multiplier_0 (
          .CLK(clk),
          .A($signed(mul_in[ IO_WIDTH*(PIXEL - i)-1 : IO_WIDTH*(PIXEL - i - 1) ])),
          .B(mul_weight),
          .P(mul_out_w[PIXEL-1-i])
        );
      end
    endgenerate
    
    



endmodule