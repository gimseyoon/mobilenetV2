//new addr_counter

`timescale 1ns / 1ps

module addr_counter (
    input clk, 
    input rst_n,
    input [2:0] state,
    input [14:0] cnt,
    input save_valid,
    input [8:0] channel_num,
// bram enable
    input in_ena,
    input in_enb,
    input weight_ena,
//input bram address
    output reg [8:0] in_addra,
    output reg [8:0] in_addrb,
//weight bram address
    output reg [14:0] weight_addra
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
                end

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
