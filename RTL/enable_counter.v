`timescale 1ns / 1ps

module enable_counter #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),        // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64),       // 15 (for 64*384 = 24576)
    parameter ADDR_W1_MEM = $clog2(384 * 9)       // 12 (for 9*384 = 3456)
)(
    input                     clk,
    input                     rst_n,
    input               [2:0] state,
    input               [3:0] bn_cnt,
    input              [14:0] glbl_cnt,
    input  [ADDR_CHANNEL-1:0] acc_cnt,
    input                     save_valid,

    // BRAM 0
    output                    ena_0,
    output                    wea_0,
    output reg                enb_0,

    // BRAM 1
    output                    ena_1,
    output                    wea_1,
    output reg                enb_1,

    // Weight BRAM
    output reg                ena_w0,
    output reg                ena_w1,
    output reg                ena_w2,
    
    // BN Parameters BRAM
    output reg                ena_bias_0,
    output reg                ena_mean_0,
    output reg                ena_std_0,
    output reg                ena_weight_0
);

////////////////////////////////////////////////////////////

    localparam IDLE         = 3'b000,
               PW_1         = 3'b001,
               PW_1_RST     = 3'b010,
               DW           = 3'b011,
               DW_RST       = 3'b100,
               PW_2         = 3'b101,
               PW_2_RST     = 3'b110,
               EXPORT       = 3'b111;

    localparam integer START   = 10'd473;   // 처음 트리거 시점
    localparam integer PERIOD  = 9'd384;    // 반복 간격
    localparam integer HIGHLEN = 2;         // high 길이(클럭 수)
    localparam integer REPEAT  = 7'd64;     // 총 반복 횟수
////////////////////////////////////////////////////////////

    reg        active;       // 반복 구간 활성화
    reg [8:0]  phase;        // 0..383 (PERIOD-1)
    reg [5:0]  rep_cnt;      // 0..63   (REPEAT-1)

///////////////////////////////////////////////////////////////////////

    assign ena_1 = ( (state==PW_1 || state==DW) && (save_valid) ) ? 1 : 0;
    assign wea_1 = ( (state==PW_1 || state==DW) && (save_valid) ) ? 1 : 0;
    assign ena_0 = ( (state==PW_2) && (save_valid) ) ? 1 : 0;
    assign wea_0 = ( (state==PW_2) && (save_valid) ) ? 1 : 0;

///////////////////////////////////////////////////////////////////////
// enable_counter

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            enb_0         <= 0;
            enb_1         <= 0;
            ena_w0        <= 0; 
            ena_w1        <= 0; 
            ena_w2        <= 0; 
            ena_bias_0    <= 0;
            ena_mean_0    <= 0;
            ena_std_0     <= 0;
            ena_weight_0  <= 0;
            active        <= 1'b0;
            phase         <= 9'd0;
            rep_cnt       <= 6'd0;
        end 
        else begin
            case (state)
                IDLE: begin
                    enb_0 <= 0; 
                    enb_1 <= 0;
                    ena_w0 <= 0; 
                    ena_w1 <= 0; 
                    ena_w2 <= 0;
                    ena_bias_0 <= 0; 
                    ena_mean_0 <= 0; 
                    ena_std_0 <= 0; 
                    ena_weight_0 <= 0;
                    active  <= 1'b0;
                    phase   <= 9'd0;
                    rep_cnt <= 6'd0;
                end
          
                PW_1: begin
                
                //bram_0, bram_W 
                    if (glbl_cnt >= 15'd24577) begin
                        enb_0  <= 0;
                        ena_w0 <= 0;
                    end
                    else begin
                        enb_0  <= 1;
                        ena_w0 <= 1;
                    end
                    
                //bram_param
                    if(bn_cnt == 4'd2) begin            
                        ena_bias_0   <= 1;
                        ena_mean_0   <= 1;  
                        ena_std_0    <= 1;
                        ena_weight_0 <= 1;                
                    end
                    else if(bn_cnt == 4'd4) begin
                        ena_bias_0   <= 0;
                        ena_mean_0   <= 0;  
                        ena_std_0    <= 0;
                        ena_weight_0 <= 0;
                    end
                end // PW_1
                
                DW: begin
                //bram_0
                    if (glbl_cnt >= 15'd5838) begin
                        enb_1  <= 0;
                        ena_w1 <= 0;
                    end
                    else begin
                        enb_1  <= 1;
                        ena_w1 <= 1;
                    end
                    
                //bram_param
                    if(glbl_cnt >= 15'd5838) begin            
                        ena_bias_0   <= 0;
                        ena_mean_0   <= 0;  
                        ena_std_0    <= 0;
                        ena_weight_0 <= 0;                
                    end
                    else begin
                        ena_bias_0   <= 1;
                        ena_mean_0   <= 1;  
                        ena_std_0    <= 1;
                        ena_weight_0 <= 1;
                    end
                end //DW

    
                PW_2: begin
                    // bram_1, bram_W2 : 원래 코드
                    if (glbl_cnt >= 15'd24577) begin
                        enb_1  <= 0;
                        ena_w2 <= 0;
                    end else begin
                        enb_1  <= 1;
                        ena_w2 <= 1;
                    end

                    if (!active) begin
                        enb_0 <= 1'b0;
                        if (glbl_cnt == START) begin
                            active  <= 1'b1;
                            phase   <= 9'd0;
                            rep_cnt <= 6'd0;
                        end
                    end else begin
                        enb_0 <= (phase < HIGHLEN);
                        
                        if (phase == (PERIOD-1)) begin
                            phase <= 9'd0;
                        if (rep_cnt == (REPEAT-1)) begin
                            // 64번 완료 후 정지
                            active  <= 1'b0;
                            enb_0   <= 1'b0;
                        end else begin
                            rep_cnt <= rep_cnt + 1'b1;
                        end
                        end else begin
                            phase <= phase + 1'b1;
                        end
                    end

                    // bram_param : 원래 코드
                    if (state==PW_2) begin
                        ena_bias_0   <= 1;
                        ena_mean_0   <= 1;  
                        ena_std_0    <= 1;
                        ena_weight_0 <= 1;
                    end else begin
                        ena_bias_0   <= 0;
                        ena_mean_0   <= 0;  
                        ena_std_0    <= 0;
                        ena_weight_0 <= 0;
                    end
                end // PW_2
    
                EXPORT: begin
                
                end
    
                default: begin
                    enb_0 <= 0; 
                    enb_1 <= 0;
                    ena_w0 <= 0; 
                    ena_w1 <= 0; 
                    ena_w2 <= 0;
                    ena_bias_0 <= 0; 
                    ena_mean_0 <= 0; 
                    ena_std_0 <= 0; 
                    ena_weight_0 <= 0;
                    active  <= 1'b0;
                    phase   <= 9'd0;
                    rep_cnt <= 6'd0;
               end
            endcase
        end
    end


endmodule