//new addr_counter

`timescale 1ns / 1ps

module addr_counter #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),        // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64),       // 15 (for 64*384 = 24576)
    parameter ADDR_W1_MEM = $clog2(384 * 9)       // 12 (for 9*384 = 3456)
)(
    input                               clk,
    input                               rst_n,
    input           [2:0]               state,
    input           [14:0]               glbl_cnt,
    input           [ADDR_CHANNEL-1:0]  acc_cnt,
    input                               save_valid,
    input                               pw_1_valid,
    input                               enb_0,
    input                               enb_1,
    output reg      [ADDR_CHANNEL-1:0] channel_num,

    // BRAM 0
    output reg      [ADDR_CHANNEL-1:0]  addra_0,
    output reg      [ADDR_CHANNEL-1:0]  addrb_0,

    // BRAM 1
    output reg      [ADDR_CHANNEL-1:0]  addra_1,
    output reg      [ADDR_CHANNEL-1:0]  addrb_1,

    // Weight BRAM
    output reg      [ADDR_WMEM-1:0]     addra_w0,
    output reg      [ADDR_W1_MEM-1:0]   addra_w1,
    output reg      [ADDR_WMEM-1:0]     addra_w2,
    
    // BN Parameter BRAMs
    output reg      [ADDR_CHANNEL-1:0]  addra_bias_0,
    output reg      [ADDR_CHANNEL-1:0]  addra_mean_0,
    output reg      [ADDR_CHANNEL-1:0]  addra_std_0,
    output reg      [ADDR_CHANNEL-1:0]  addra_weight_0
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

    reg [3:0] cnt_9;
    
////////////////////////////////////////////////////////////
// addr_counter


    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            cnt_9           <= 0;
            channel_num     <= 0;
            addra_0         <= 0;
            addrb_0         <= 0;
            addra_1         <= 0;
            addrb_1         <= 0;
            addra_w0        <= 0;
            addra_w1        <= 0;
            addra_w2        <= 0;
            addra_bias_0    <= 0;
            addra_mean_0    <= 0;
            addra_std_0     <= 0;
            addra_weight_0  <= 0;
        end 
        else begin
            case (state)
                IDLE: begin
                    cnt_9           <= 0;
                    channel_num     <= 0;
                    addra_0         <= 0;
                    addrb_0         <= 0;
                    addra_1         <= 0;
                    addrb_1         <= 0;
                    addra_w0        <= 0;
                    addra_w1        <= 0;
                    addra_w2        <= 0;
                    addra_bias_0    <= 0;
                    addra_mean_0    <= 0;
                    addra_std_0     <= 0;
                    addra_weight_0  <= 0;
                end

                PW_1: begin
                    if(enb_0) begin
                        if (glbl_cnt < 15'd24576) begin
                            addrb_0 <= (addrb_0 >= 9'd63) ? 0 : addrb_0 + 1;
                            addra_w0 <= addra_w0 + 1;
                        end else begin
                            addrb_0  <= 0;
                            addra_w0 <= 0;
                        end
                    end
                    
                    if(acc_cnt == 9'd63) begin
                        addra_1 <= addra_1 + 1;
                    end
                // bram_param
                    if(pw_1_valid) begin
                        addra_bias_0 <= addra_bias_0 + 1; 
                        addra_mean_0 <= addra_mean_0 + 1; 
                        addra_std_0 <= addra_std_0 + 1; 
                        addra_weight_0 <= addra_weight_0 + 1; 
                    end
                    
                end //PW_1
                
                PW_1_BN_RELU: begin

                end
                DW: begin
                    if(enb_1) begin
                        if (glbl_cnt < 15'd3456) begin
                            if (cnt_9 == 4'd8) begin
                                cnt_9 <= 0;
                                addrb_1 <= addrb_1 + 1;
                            end 
                            else begin
                                cnt_9 <= cnt_9 + 1;
                            end
                            addra_w1 <= addra_w1 + 1;
                        end 
                        else begin
                            addrb_1  <= 0;
                            addra_w1 <= 0;
                        end
                    end
                end
                // Other states: modify in future
                default: begin
                    cnt_9           <= 0;
                    channel_num     <= 0;
                    addra_0         <= 0;
                    addrb_0         <= 0;
                    addra_1         <= 0;
                    addrb_1         <= 0;
                    addra_w0        <= 0;
                    addra_bias_0    <= 0;
                    addra_mean_0    <= 0;
                    addra_std_0     <= 0;
                    addra_weight_0  <= 0;
                end
            endcase
        end
    end

endmodule