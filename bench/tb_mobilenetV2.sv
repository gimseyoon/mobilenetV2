`timescale 1ns / 1ps

module tb_mobilenetV2 #(

    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),        // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64),       // 15 (for 64*384 = 24576)
    parameter ADDR_W1_MEM = $clog2(384 * 9)       // 12 (for 9*384 = 3456)
);

    reg clk;
    reg rst_n;
    reg start;
    wire [13:0] layer_8_result;


    
    
/////////////////////////////////////////////////////////////////
    
    // DUT instantiation
    mobilenetV2 #(
        .IO_WIDTH(IO_WIDTH),
        .ROW(ROW),
        .COLUMN(COLUMN),
        .PIXEL(PIXEL),
        .W_WIDTH(W_WIDTH),
        .ADDR_CHANNEL(ADDR_CHANNEL),
        .ADDR_WMEM(ADDR_WMEM)
    ) dut (
        .clk(clk),
        .rst_n(rst_n),
        .start(start),
        .layer_8_result(layer_8_result)
    );

    // Clock generation
    initial clk = 0;
    always #1.6666667 clk = ~clk;

    // rst_n
    initial begin
        rst_n = 0;
        #25 rst_n = 1;
    end

    // Start
    initial begin
        start = 0;
        #400 start = 1;
        #25 start = 0;
    end


        
endmodule