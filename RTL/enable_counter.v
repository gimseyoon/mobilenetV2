`timescale 1ns / 1ps

module enable_counter #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),         // 8 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64)       // 15 (for 64*384 = 24576)
)(
    input clk,
    input rst_n,
    input [2:0] state,
    input [ADDR_WMEM-1 : 0] cnt,
// input bram eanble
    output reg in_ena,
    output reg in_wea,
    output reg in_enb,
// weight bram enable
    output reg weight_ena
);

////////////////////////////////////////////////////////////

    localparam IDLE         = 3'b000,
               PW_1         = 3'b001,
               PW_1_BN_RELU = 3'b010,
               DW           = 3'b011,
               DW_BN_RELU   = 3'b100,
               PW_2         = 3'b101,
               PW_2_BN      = 3'b110,
               SK           = 3'b111;

///////////////////////////////////////////////////////////////////////
// enable_counter

always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            in_ena      <= 0;
            in_wea      <= 0;
            in_enb      <= 0;
            weight_ena  <= 0;
        end 
        else begin
            // 기본값 초기화
            in_ena      <= 0;
            in_wea      <= 0;
            in_enb      <= 0;
            weight_ena  <= 0;

            case (state)
                PW_1: begin
                    if (cnt >= 15'd24577) begin
                        in_enb     <= 0;
                        weight_ena <= 0;
                    end
                    else begin
                        in_enb     <= 1;
                        weight_ena <= 1;
                    end
                end // begin

                PW_1_BN_RELU: begin
                end

                DW: begin
                    in_enb     <= 1;
                    weight_ena <= 1;
                end

                // 나머지 상태는 기본값 유지 (0)
            endcase
        end //else
    end //always

endmodule