`timescale 1ns / 1ps

module SK #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter INPUT_CHANNEL = 64,
    parameter ADDR_PARAM = 10,
    parameter ADDR_CHANNEL  = $clog2(384),        // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64),       // 15 (for 64*384 = 24576)
    parameter ADDR_W1_MEM = $clog2(384 * 9)       // 12 (for 9*384 = 3456)
)(
    input                                   clk,
    input                                   rst_n,
    input                                   skip_valid,
    input       signed [IO_WIDTH*PIXEL-1:0] in_1,
    input       signed [IO_WIDTH*PIXEL-1:0] in_2,
    output      signed [IO_WIDTH*PIXEL-1:0] result,
    output reg                              result_save_valid
    );
    
reg                         skip_valid_d;
reg  signed [IO_WIDTH*PIXEL-1:0]   in_1_q;
reg  signed [IO_WIDTH*PIXEL-1:0]   in_2_q;
wire signed [IO_WIDTH-1:0]         in_1_reg    [0:PIXEL-1];
wire signed [IO_WIDTH-1:0]         in_2_reg    [0:PIXEL-1];
wire signed [IO_WIDTH-1:0]         sum_reg     [0:PIXEL-1];
reg  signed [IO_WIDTH-1:0]         result_reg  [0:PIXEL-1];

/////////////////////////////////////////////////////// 
// Input bus pipeline ( for short net delay )
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        in_1_q <= 0;
        in_2_q <= 0;
        skip_valid_d <= 0;
    end else begin
        in_1_q <= in_1;   
        in_2_q <= in_2;
        skip_valid_d <= skip_valid;
    end
end

genvar j;
generate
    for (j = 0; j < PIXEL; j = j + 1) begin : SLICE
        assign in_1_reg[j] = in_1_q[IO_WIDTH*(j+1)-1 : IO_WIDTH*j];
        assign in_2_reg[j] = in_2_q[IO_WIDTH*(j+1)-1 : IO_WIDTH*j];
        assign sum_reg[j] = in_1_reg[j] + in_2_reg[j]; // 고정소수점 동일 포맷 가정
    end
endgenerate
    

integer k;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        result_save_valid <= 0;
        for (k = 0; k < PIXEL; k = k + 1) begin
            result_reg[k] <= 0;
        end
    end else begin
        if (skip_valid_d) begin
            result_save_valid <= 1;
            for (k = 0; k < PIXEL; k = k + 1) begin
                result_reg[k] <= sum_reg[k];
            end
        end
        else begin
            result_save_valid <= 0;
        end
    end
end
    
/////////////////////////////////////////////////////// 
// Pack result_reg[] -> result (bus mapping)
///////////////////////////////////////////////////////
genvar m;
generate
    for (m = 0; m < PIXEL; m = m + 1) begin : PACK_RESULT
        assign result[IO_WIDTH*(m+1)-1 : IO_WIDTH*m] = result_reg[m];
    end
endgenerate
    
    
    
endmodule
