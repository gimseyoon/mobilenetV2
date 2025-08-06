//new addr_counter

`timescale 1ns / 1ps

module addr_counter #(
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
    input [ADDR_WMEM-1 :0] cnt,
    input [ADDR_CHANNEL -1 :0]  acc_cnt,
    input save_valid,
    input [ADDR_CHANNEL-1 : 0] channel_num,
// bram enable
    input in_ena,
    input in_enb,
    input weight_ena,
    
//input bram address
    output reg [ADDR_CHANNEL-1 : 0] in_addra,
    output reg [ADDR_CHANNEL-1 : 0] in_addrb,
//weight bram address
    output reg [ADDR_WMEM-1 : 0] weight_addra
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

////////////////////////////////////////////////////////////
// addr_counter


    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            in_addra     <= 0;
            in_addrb     <= 0;
            weight_addra <= 0;
        end else begin
            case (state)
                IDLE: begin
                    in_addra     <= 0;
                    in_addrb     <= 0;
                    weight_addra <= 0;
                end

                PW_1: begin
                    if(in_enb) begin
                        if (cnt < 15'd24576) begin
                            in_addrb     <= (in_addrb >= 9'd63) ? 0 : in_addrb + 1;
                            weight_addra <= weight_addra + 1;
                        end else begin
                            in_addrb     <= 0;
                            weight_addra <= 0;
                        end
                    end
                    
                    
                end //PW_1

                // Other states: modify in future
                default: begin
                    in_addra     <= 0;
                    in_addrb     <= 0;
                    weight_addra <= 0;
                end
            endcase
        end
    end

endmodule