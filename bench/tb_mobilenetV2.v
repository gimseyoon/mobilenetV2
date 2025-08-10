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
    reg [31:0] test_value_1 = 32'hFFFEF2F8; // 1111_1111_1111_1110_1111_0010_1111_1000 ->  110_1111_0010_1111_100
    reg [17:0] test_value_2 = 18'b110_1110_0101_1111_000;
    reg [17:0] test_value_3 = 18'b110_1111_0010_1111_100;
    reg [19:0] test_value_4 = 20'hEF2F8; // 1110_1111_0010_1111_1000
    
    
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
        .start(start)
    );

    // Clock generation
    initial clk = 0;
    always #12.5 clk = ~clk;

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