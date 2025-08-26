`timescale 1ns / 1ps
module reset_sync #(
    parameter integer STAGES = 2 
)(
    input  wire clk,
    input  wire rst,     
    output wire rst_n_sync, 
    output wire rst_sync     
);

    initial begin
        if (STAGES < 2) $error("reset_sync: STAGES must be >= 2");
    end

    (* ASYNC_REG = "TRUE" *) reg [STAGES-1:0] shreg;


    always @(posedge clk or posedge rst) begin
        if (rst) begin

            shreg <= {STAGES{1'b0}};
        end else begin

            shreg <= {shreg[STAGES-2:0], 1'b1};
        end
    end

    assign rst_n_sync = shreg[STAGES-1];      
    assign rst_sync   = ~rst_n_sync;         
endmodule