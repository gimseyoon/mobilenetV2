
`timescale 1ns / 1ps

module addr_counter #(
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
    input                               clk,
    input                               rst_n,
    input           [2:0]               state,
    input           [3:0]               bn_cnt,
    input           [14:0]              glbl_cnt,
    input           [ADDR_CHANNEL-1:0]  acc_cnt,
    input                       [12:0]  dw_cnt, // 0 ~ 5760( = (14+1) x 384 )
    input                               save_valid,
    input                               skip_valid,
    input                               enb_0,
    input                               enb_1,

    // BRAM 0
    output reg     [INPUT_CHANNEL-1:0]  addra_0,
    output reg     [INPUT_CHANNEL-1:0]  addrb_0,

    // BRAM 1
    output reg      [ADDR_CHANNEL-1:0]  addra_1,
    output reg      [ADDR_CHANNEL-1:0]  addrb_1,

    // Weight BRAM
    output reg      [ADDR_WMEM-1:0]     addra_w0,
    output reg      [ADDR_W1_MEM-1:0]   addra_w1,
    output reg      [ADDR_WMEM-1:0]     addra_w2,
    
    // BN Parameter BRAMs
    output reg        [ADDR_PARAM-1:0]  addra_bias_0,
    output reg        [ADDR_PARAM-1:0]  addra_mean_0,
    output reg        [ADDR_PARAM-1:0]  addra_std_0,
    output reg        [ADDR_PARAM-1:0]  addra_weight_0
);
////////////////////////////////////////////////////////////

    localparam IDLE         = 3'b000,
               PW_1         = 3'b001,
               PW_1_RST     = 3'b010,
               DW           = 3'b011,
               DW_RST       = 3'b100,
               PW_2         = 3'b101,
               PW_2_RST     = 3'b110,
               EXPORT       = 3'b111,
               
               DW_OFFSET    = 10'd384,
               PW_2_OFFSET  = 10'd768;
                
////////////////////////////////////////////////////////////
    reg enb_0_q;
    reg [3:0] cnt;
    reg       mean_run,  std_run,  weight_run,  bias_run;
    reg [3:0] mean_phase,std_phase,weight_phase,bias_phase;
////////////////////////////////////////////////////////////
// addr_counter

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            enb_0_q          <= 0;
            cnt             <= 0;
            addra_0         <= 0;
            addrb_0         <= 0;
            addra_1         <= 0;
            addrb_1         <= 0;
            addra_w0        <= 0;
            addra_w1        <= 0;
            addra_w2        <= 0;
            addra_bias_0    <= 10'd1023;
            addra_mean_0    <= 10'd1023;
            addra_std_0     <= 10'd1023;
            addra_weight_0  <= 10'd1023;
            mean_run   <= 1'b0; std_run   <= 1'b0; weight_run   <= 1'b0; bias_run   <= 1'b0;
            mean_phase <= 4'd0; std_phase <= 4'd0; weight_phase <= 4'd0; bias_phase <= 4'd0;
        end 
        else begin
            enb_0_q <= enb_0;
            case (state)
                IDLE: begin
                    enb_0_q          <= 0;
                    cnt             <= 0;
                    addra_0         <= 0;
                    addrb_0         <= 0;
                    addra_1         <= 0;
                    addrb_1         <= 0;
                    addra_w0        <= 0;
                    addra_w1        <= 0;
                    addra_w2        <= 0;
                    addra_bias_0    <= 10'd1023;
                    addra_mean_0    <= 10'd1023;
                    addra_std_0     <= 10'd1023;
                    addra_weight_0  <= 10'd1023;
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
                    
                    if(save_valid) begin
                        addra_1 <= addra_1 + 1;
                    end
                // bram_param
                    if(bn_cnt == 4'd2) begin
                        addra_bias_0 <= addra_bias_0 + 1; 
                        addra_mean_0 <= addra_mean_0 + 1; 
                        addra_std_0 <= addra_std_0 + 1; 
                        addra_weight_0 <= addra_weight_0 + 1; 
                    end
                    
                end //PW_1
                
                PW_1_RST: begin
                    addra_1 <= 0;
                    addra_bias_0    <= DW_OFFSET;
                    addra_mean_0    <= DW_OFFSET;
                    addra_std_0     <= DW_OFFSET;
                    addra_weight_0  <= DW_OFFSET;
                    mean_run   <= 1'b0; std_run   <= 1'b0; weight_run   <= 1'b0; bias_run   <= 1'b0;
                    mean_phase <= 4'd0; std_phase <= 4'd0; weight_phase <= 4'd0; bias_phase <= 4'd0;
                end
                
                DW: begin
                    if(enb_1) begin
                        if (glbl_cnt < 15'd5838) begin
                            if((cnt <= 4'd7) || (cnt == 4'd14)) begin
                                if(addra_w1 == 3455) begin
                                    addra_w1 <= 0;
                                end
                                else begin
                                    addra_w1 <= addra_w1 + 1;
                                end

                            end
                            
                            if (cnt == 4'd14) begin
                                cnt <= 0;
                                if(addrb_1 == 383) begin
                                    addrb_1 <= 0;
                                end
                                else begin
                                    addrb_1 <= addrb_1 + 1;
                                end

                            end 
                            else begin
                                cnt <= cnt + 1;
                            end

                        end 
                        else begin
                            addrb_1  <= 0;
                            addra_w1 <= 0;
                        end
                    end // if(enb_1)
                    
                    if(save_valid) begin
                        if(addra_1 == 384) addra_1 <= 0;
                        else               addra_1 <= addra_1 + 1;
                    end

                // ===== mean: base=36 이후 15주기 =====
                    if (glbl_cnt == 13'd35) begin
                        addra_mean_0 <= addra_mean_0 + 1'b1;   // 베이스에서 1회 증가
                        mean_run     <= 1'b1;                  // 런 시작
                        mean_phase   <= 4'd0;                  // 주기 카운터 0으로 정렬
                    end 
                    else if (mean_run) begin
                        if (mean_phase == 4'd14) begin         // 15클럭마다
                            addra_mean_0 <= addra_mean_0 + 1'b1;
                            mean_phase   <= 4'd0;
                        end 
                        else begin
                        mean_phase <= mean_phase + 1'b1;
                        end
                    end
              
                    // ===== std: base=47 이후 15주기 =====
                    if (glbl_cnt == 13'd46) begin
                        addra_std_0  <= addra_std_0 + 1'b1;
                        std_run      <= 1'b1;
                        std_phase    <= 4'd0;
                    end 
                    else if (std_run) begin
                        if (std_phase == 4'd14) begin
                            addra_std_0 <= addra_std_0 + 1'b1;
                            std_phase   <= 4'd0;
                        end else begin
                            std_phase <= std_phase + 1'b1;
                        end
                    end
              
                    // ===== weight: base=75 이후 15주기 =====
                    if (glbl_cnt == 13'd74) begin
                        addra_weight_0 <= addra_weight_0 + 1'b1;
                        weight_run     <= 1'b1;
                        weight_phase   <= 4'd0;
                    end else if (weight_run) begin
                        if (weight_phase == 4'd14) begin
                            addra_weight_0 <= addra_weight_0 + 1'b1;
                            weight_phase   <= 4'd0;
                        end else begin
                            weight_phase <= weight_phase + 1'b1;
                        end
                    end
              
                    // ===== bias: base=83 이후 15주기 =====
                    if (glbl_cnt == 13'd82) begin
                        addra_bias_0 <= addra_bias_0 + 1'b1;
                        bias_run     <= 1'b1;
                        bias_phase   <= 4'd0;
                    end else if (bias_run) begin
                        if (bias_phase == 4'd14) begin
                            addra_bias_0 <= addra_bias_0 + 1'b1;
                            bias_phase   <= 4'd0;
                        end else begin
                            bias_phase <= bias_phase + 1'b1;
                        end
                    end       
                    
                end // DW
                
                DW_RST: begin
                    addra_1 <= 0;
                    addra_bias_0    <= PW_2_OFFSET;
                    addra_mean_0    <= PW_2_OFFSET;
                    addra_std_0     <= PW_2_OFFSET;
                    addra_weight_0  <= PW_2_OFFSET;
                    mean_run   <= 1'b0; std_run   <= 1'b0; weight_run   <= 1'b0; bias_run   <= 1'b0;
                    mean_phase <= 4'd0; std_phase <= 4'd0; weight_phase <= 4'd0; bias_phase <= 4'd0;
                end
                
                PW_2: begin
                    if(enb_1) begin
                        if (glbl_cnt < 15'd24576) begin
                            addrb_1 <= (addrb_1 >= 9'd383) ? 0 : addrb_1 + 1;
                            addra_w2 <= addra_w2 + 1;
                        end 
                        else begin
                            addrb_1  <= 0;
                            addra_w2 <= 0;
                        end
                    end
                    
                    // bram_param
                    if(skip_valid) begin
                        addra_1 <= addra_1 + 1;
                        addra_bias_0 <= addra_bias_0 + 1; 
                        addra_mean_0 <= addra_mean_0 + 1; 
                        addra_std_0 <= addra_std_0 + 1; 
                        addra_weight_0 <= addra_weight_0 + 1; 
                    end

                    if (enb_0_q && !enb_0) begin
                        if (addrb_0 == 6'd63)  // 필요시 9'd63 등으로 바꿔도 OK
                            addrb_0 <= 0;
                        else
                            addrb_0 <= addrb_0 + 1'b1;
                    end
                        
                        
                end // PW_2
                                
                EXPORT: begin
                
                end // SK
                
                default: begin
                    cnt           <= 0;
                    addra_0         <= 0;
                    addrb_0         <= 0;
                    addra_1         <= 0;
                    addrb_1         <= 0;
                    addra_w0        <= 0;
                    addra_w1        <= 0;
                    addra_w2        <= 0;
                    addra_bias_0    <= 10'd1023;
                    addra_mean_0    <= 10'd1023;
                    addra_std_0     <= 10'd1023;
                    addra_weight_0  <= 10'd1023;
                end // Default
                
            endcase
        end
    end

endmodule