`timescale 1ns / 1ps

module tb_mobilenetV2#(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW*COLUMN, // 14 * 14 = 196
    parameter W_WIDTH = 17
);

    reg clk;
    reg rst_n;
    reg start;

    
    
/////////////////////////////////////////////////////////////////
    
    // DUT instantiation
    mobilenetV2 #(
        .IO_WIDTH(IO_WIDTH),
        .ROW(ROW),
        .COLUMN(COLUMN),
        .PIXEL(PIXEL),
        .W_WIDTH(W_WIDTH)
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
        rst_n = 1;
        #300 rst_n = 0;
        #25 rst_n = 1;
    end

    // Start
    initial begin
        start = 0;
        #400 start = 1;
        #25 start = 0;
    end


        
endmodule