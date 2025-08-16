`timescale 1ns / 1ps

module addr_counter #(
    parameter IO_WIDTH       = 18,
    parameter ROW            = 14,
    parameter COLUMN         = 14,
    parameter PIXEL          = ROW * COLUMN,              // 196
    parameter W_WIDTH        = 17,
    parameter INPUT_CHANNEL  = 64,
    parameter ADDR_PARAM     = 10,
    parameter ADDR_CHANNEL   = $clog2(384),               // 9
    parameter ADDR_WMEM      = $clog2(384 * 64),          // 15
    parameter ADDR_W1_MEM    = $clog2(384 * 9)            // 12
)(
    input                           clk,
    input                           rst_n,
    input         [2:0]             state,
    input         [14:0]            glbl_cnt,
    input                           save_valid,
    input                           skip_valid,
    input                           enb_0,
    input                           enb_1,

    output reg                     pw_1_read_done,
    output reg                     dw_read_done,
    output reg                     pw_2_read_done,
    // BRAM 0
    output reg    [INPUT_CHANNEL-1:0] addra_0,
    output reg    [INPUT_CHANNEL-1:0] addrb_0,

    // BRAM 1
    output reg    [ADDR_CHANNEL-1:0]  addra_1,
    output reg    [ADDR_CHANNEL-1:0]  addrb_1,

    // Weights
    output reg    [ADDR_WMEM-1:0]     addra_w0,
    output reg    [ADDR_W1_MEM-1:0]   addra_w1,
    output reg    [ADDR_WMEM-1:0]     addra_w2,

    // BN params
    output reg    [ADDR_PARAM-1:0]    addra_bias_0,
    output reg    [ADDR_PARAM-1:0]    addra_mean_0,
    output reg    [ADDR_PARAM-1:0]    addra_std_0,
    output reg    [ADDR_PARAM-1:0]    addra_weight_0
);

///////////////////////////////////////////////////////
// States & offsets
///////////////////////////////////////////////////////
localparam IDLE       = 3'b000,
           PW_1       = 3'b001,
           PW_1_RST   = 3'b010,
           DW         = 3'b011,
           DW_RST     = 3'b100,
           PW_2       = 3'b101,
           PW_2_RST   = 3'b110,
           EXPORT     = 3'b111;

localparam [ADDR_PARAM-1:0] DW_OFFSET   = 10'd384;
localparam [ADDR_PARAM-1:0] PW_2_OFFSET = 10'd768;

///////////////////////////////////////////////////////
// Registers
///////////////////////////////////////////////////////
reg         enb_0_q;
reg  [4:0]  dw_cnt;
reg         mean_run, std_run, weight_run, bias_run;
reg  [5:0]  mean_phase, std_phase, weight_phase, bias_phase;

///////////////////////////////////////////////////////
// Sequential: address generation
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        enb_0_q        <= 1'b0;
        dw_cnt            <= 5'd0;

        addra_0        <= {INPUT_CHANNEL{1'b0}};
        addrb_0        <= {INPUT_CHANNEL{1'b0}};
        addra_1        <= {ADDR_CHANNEL{1'b0}};
        addrb_1        <= {ADDR_CHANNEL{1'b0}};
        addra_w0       <= {ADDR_WMEM{1'b0}};
        addra_w1       <= {ADDR_W1_MEM{1'b0}};
        addra_w2       <= {ADDR_WMEM{1'b0}};

        addra_bias_0   <= {ADDR_PARAM{1'b0}};  // 1023
        addra_mean_0   <= {ADDR_PARAM{1'b0}};
        addra_std_0    <= {ADDR_PARAM{1'b0}};
        addra_weight_0 <= {ADDR_PARAM{1'b0}};

        mean_run       <= 1'b0;  std_run   <= 1'b0;  weight_run <= 1'b0;  bias_run <= 1'b0;
        mean_phase     <= 6'd0;  std_phase <= 6'd0;  weight_phase<= 6'd0; bias_phase<= 6'd0;
    end else begin
        enb_0_q <= enb_0;

        case (state)
            ///////////////////////////////////////////////////////
            // IDLE
            ///////////////////////////////////////////////////////
            IDLE: begin
                enb_0_q         <= 1'b0;
                dw_cnt          <= 5'd0;
                pw_1_read_done  <= 0;
                dw_read_done    <= 0;
                pw_2_read_done  <= 0;
                addra_0         <= {INPUT_CHANNEL{1'b0}};
                addrb_0         <= {INPUT_CHANNEL{1'b0}};
                addra_1         <= {ADDR_CHANNEL{1'b0}};
                addrb_1         <= {ADDR_CHANNEL{1'b0}};
                addra_w0        <= {ADDR_WMEM{1'b0}};
                addra_w1        <= {ADDR_W1_MEM{1'b0}};
                addra_w2        <= {ADDR_WMEM{1'b0}};

                addra_bias_0    <= {ADDR_PARAM{1'b0}};
                addra_mean_0    <= {ADDR_PARAM{1'b0}};
                addra_std_0     <= {ADDR_PARAM{1'b0}};
                addra_weight_0  <= {ADDR_PARAM{1'b0}};
            end

            ///////////////////////////////////////////////////////
            // PW_1
            ///////////////////////////////////////////////////////
            PW_1: begin
                // BRAM_A, BRAM_W_0
                if (enb_0) begin
                    addrb_0  <= (addrb_0 >= 9'd63) ? {INPUT_CHANNEL{1'b0}} : addrb_0 + 1'b1;
                    addra_w0 <= addra_w0 + 1'b1;
                end
                else begin
                    addrb_0  <= {INPUT_CHANNEL{1'b0}};
                    addra_w0 <= {ADDR_WMEM{1'b0}};
                end

                if(addra_w0 == 15'd24576) begin
                    pw_1_read_done <= 1;
                end
                
                if (save_valid)
                    addra_1 <= addra_1 + 1'b1;
            
                // -------------------------------------------------
                // BN params (PW_1): base 트리거 후 매 64클럭 증가
                //   - mean  : glbl_cnt == 105 -> 다음 clk에 1증가, 이후 매 64클럭
                //   - std   : glbl_cnt == 116 -> 다음 clk에 1증가, 이후 매 64클럭
                //   - weight: glbl_cnt == 144 -> 다음 clk에 1증가, 이후 매 64클럭
                //   - bias  : glbl_cnt == 152 -> 다음 clk에 1증가, 이후 매 64클럭
                // -------------------------------------------------
                // mean
                if (glbl_cnt == 15'd106) begin
                    mean_run   <= 1'b1;
                    mean_phase <= 6'd63;           
                end else if (mean_run) begin
                    if (mean_phase == 6'd63) begin
                        addra_mean_0 <= addra_mean_0 + 1'b1;
                        mean_phase   <= 6'd0;
                    end else begin
                        mean_phase   <= mean_phase + 1'b1;
                    end
                end
            
                // std
                if (glbl_cnt == 15'd117) begin
                    std_run    <= 1'b1;
                    std_phase  <= 6'd63;
                end else if (std_run) begin
                    if (std_phase == 6'd63) begin
                        addra_std_0 <= addra_std_0 + 1'b1;
                        std_phase   <= 6'd0;
                    end else begin
                        std_phase   <= std_phase + 1'b1;
                    end
                end
            
                // weight
                if (glbl_cnt == 15'd145) begin
                    weight_run   <= 1'b1;
                    weight_phase <= 6'd63;
                end else if (weight_run) begin
                    if (weight_phase == 6'd63) begin
                        addra_weight_0 <= addra_weight_0 + 1'b1;
                        weight_phase   <= 6'd0;
                    end else begin
                        weight_phase   <= weight_phase + 1'b1;
                    end
                end
            
                // bias
                if (glbl_cnt == 15'd153) begin
                    bias_run   <= 1'b1;
                    bias_phase <= 6'd63;
                end else if (bias_run) begin
                    if (bias_phase == 6'd63) begin
                        addra_bias_0 <= addra_bias_0 + 1'b1;
                        bias_phase   <= 6'd0;
                    end else begin
                        bias_phase   <= bias_phase + 1'b1;
                    end
                end
            end // PW_1

            ///////////////////////////////////////////////////////
            // PW_1_RST
            ///////////////////////////////////////////////////////
            PW_1_RST: begin
                addra_1        <= {ADDR_CHANNEL{1'b0}};
                addra_bias_0   <= DW_OFFSET;
                addra_mean_0   <= DW_OFFSET;
                addra_std_0    <= DW_OFFSET;
                addra_weight_0 <= DW_OFFSET;

                mean_run<=1'b0; std_run<=1'b0; weight_run<=1'b0; bias_run<=1'b0;
                mean_phase<=6'd0; std_phase<=6'd0; weight_phase<=6'd0; bias_phase<=6'd0;
            end

            ///////////////////////////////////////////////////////
            // DW
            ///////////////////////////////////////////////////////
            DW: begin

                // BRAM_B, BRAM_W1
                if (enb_1) begin
                        if ((dw_cnt <= 5'd7) || (dw_cnt == 5'd28)) begin
                            addra_w1 <= addra_w1 + 1'b1;
                        end
                        if (dw_cnt == 5'd28) begin
                            dw_cnt     <= 4'd0;
                            addrb_1 <= (addrb_1 == 9'd383) ? {ADDR_CHANNEL{1'b0}} : addrb_1 + 1'b1;
                        end else begin
                            dw_cnt <= dw_cnt + 1'b1;
                        end
                end
                else begin
                    addra_w1 <= 0;
                    addrb_1 <= 0;
                end
                
                if(addra_w1 == 12'd3455) begin
                    dw_read_done <= 1;
                end
                
                if (save_valid)
                    addra_1 <= addra_1 + 1'b1;  

                // mean: base=35+1, every 15
                if (glbl_cnt == 15'd49) begin
                    addra_mean_0 <= addra_mean_0 + 1'b1; mean_run<=1'b1; mean_phase<=6'd0;
                end else if (mean_run) begin
                    if (mean_phase == 6'd28) begin
                        addra_mean_0 <= addra_mean_0 + 1'b1; mean_phase <= 6'd0;
                    end else begin
                        mean_phase <= mean_phase + 6'b1;
                    end
                end

                // std: base=46, every 15
                if (glbl_cnt == 15'd60) begin
                    addra_std_0 <= addra_std_0 + 1'b1; std_run<=1'b1; std_phase<=4'd0;
                end else if (std_run) begin
                    if (std_phase == 6'd28) begin
                        addra_std_0 <= addra_std_0 + 1'b1; std_phase <= 4'd0;
                    end else begin
                        std_phase <= std_phase + 1'b1;
                    end
                end

                // weight: base=74, every 15
                if (glbl_cnt == 15'd88) begin
                    addra_weight_0 <= addra_weight_0 + 1'b1; weight_run<=1'b1; weight_phase<=4'd0;
                end else if (weight_run) begin
                    if (weight_phase == 6'd28) begin
                        addra_weight_0 <= addra_weight_0 + 1'b1; weight_phase <= 4'd0;
                    end else begin
                        weight_phase <= weight_phase + 1'b1;
                    end
                end

                // bias: base=82, every 15
                if (glbl_cnt == 15'd96) begin
                    addra_bias_0 <= addra_bias_0 + 1'b1; bias_run<=1'b1; bias_phase<=4'd0;
                end else if (bias_run) begin
                    if (bias_phase == 6'd28) begin
                        addra_bias_0 <= addra_bias_0 + 1'b1; bias_phase <= 4'd0;
                    end else begin
                        bias_phase <= bias_phase + 1'b1;
                    end
                end
            end

            ///////////////////////////////////////////////////////
            // DW_RST
            ///////////////////////////////////////////////////////
            DW_RST: begin
                addra_1        <= {ADDR_CHANNEL{1'b0}};
                addra_bias_0   <= PW_2_OFFSET;
                addra_mean_0   <= PW_2_OFFSET;
                addra_std_0    <= PW_2_OFFSET;
                addra_weight_0 <= PW_2_OFFSET;

                mean_run<=1'b0; std_run<=1'b0; weight_run<=1'b0; bias_run<=1'b0;
                mean_phase<=6'd0; std_phase<=6'd0; weight_phase<=6'd0; bias_phase<=6'd0;
            end

            ///////////////////////////////////////////////////////
            // PW_2
            ///////////////////////////////////////////////////////
            PW_2: begin
                // BRAM_B, BRAM_W2
                if (enb_1) begin
                    addrb_1  <= (addrb_1 >= 9'd383) ? {ADDR_CHANNEL{1'b0}} : addrb_1 + 1'b1;
                    addra_w2 <= addra_w2 + 1'b1;
                end 
                else begin
                    addrb_1  <= {ADDR_CHANNEL{1'b0}};
                    addra_w2 <= {ADDR_WMEM{1'b0}};
                end
                    
                if(addra_w2 == 15'd24575) begin
                    pw_2_read_done <= 1;
                end
               

                                                                
                // PW2 uses skip_valid
                if (skip_valid) begin
                    addra_1        <= addra_1 + 1'b1;
                    addra_bias_0   <= addra_bias_0   + 1'b1;
                    addra_mean_0   <= addra_mean_0   + 1'b1;
                    addra_std_0    <= addra_std_0    + 1'b1;
                    addra_weight_0 <= addra_weight_0 + 1'b1;
                end

                // addra_0++ on enb_0 falling edge (64 times total)
                if (enb_0_q && !enb_0)
                    addra_0 <= (addra_0 == 6'd63) ? {INPUT_CHANNEL{1'b0}} : addra_0 + 1'b1;
            end

            ///////////////////////////////////////////////////////
            // EXPORT (reserved)
            ///////////////////////////////////////////////////////
            EXPORT: begin
                // no-op
            end

            ///////////////////////////////////////////////////////
            // default
            ///////////////////////////////////////////////////////
            default: begin
                pw_1_read_done <= 0;
                dw_read_done <= 0;
                pw_2_read_done <= 0;
                dw_cnt         <= 5'd0;
                addra_0        <= {INPUT_CHANNEL{1'b0}};
                addrb_0        <= {INPUT_CHANNEL{1'b0}};
                addra_1        <= {ADDR_CHANNEL{1'b0}};
                addrb_1        <= {ADDR_CHANNEL{1'b0}};
                addra_w0       <= {ADDR_WMEM{1'b0}};
                addra_w1       <= {ADDR_W1_MEM{1'b0}};
                addra_w2       <= {ADDR_WMEM{1'b0}};
                addra_bias_0   <= {ADDR_PARAM{1'b0}};
                addra_mean_0   <= {ADDR_PARAM{1'b0}};
                addra_std_0    <= {ADDR_PARAM{1'b0}};
                addra_weight_0 <= {ADDR_PARAM{1'b0}};
            end
        endcase
    end
end

endmodule
