`timescale 1ns / 1ps

module addr_counter #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter INPUT_CHANNEL = 64,
    parameter ADDR_PARAM = 12,
    parameter ADDR_IN = $clog2(64),
    parameter ADDR_CHANNEL  = $clog2(384),        // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64),       // 15 (for 64*384 = 24576)
    parameter ADDR_W1_MEM = $clog2(384 * 9)       // 12 (for 9*384 = 3456)
)(
    input                           clk,
    input                           rst_n,
    input         [2:0]             state,
    input         [1:0]             layer_state,
    input                           save_valid,
    input                           skip_valid,
    input                           result_save_valid,
    input                           enb_0,
    input                           enb_1,

    output reg                     pw_1_read_done,
    output reg                     dw_read_done,
    output reg                     pw_2_read_done,
    // BRAM 0
    output reg    [ADDR_IN-1:0]     addra_0,
    output reg    [ADDR_IN-1:0]     addrb_0,

    // BRAM 1
    output reg    [ADDR_CHANNEL-1:0]  addra_1,
    output reg    [ADDR_CHANNEL-1:0]  addrb_1,

    // Weights
    output reg    [ADDR_WMEM-1:0]     addra_w0,
    output reg    [ADDR_W1_MEM-1:0]   addra_w1,
    output reg    [ADDR_WMEM-1:0]     addra_w2,

    // BN params
    output reg    [ADDR_PARAM-1:0]    addra_biassubam_0,
    output reg    [ADDR_PARAM-1:0]    addra_wdivstd_0
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

localparam READY    = 2'b00,
           LAYER_8  = 2'b01,
           LAYER_9  = 2'b10,
           LAYER_10 = 2'b11;
           
localparam [ADDR_PARAM-1:0] LAYER_8_PW_1_OFFSET  = 12'd0;
localparam [ADDR_PARAM-1:0] LAYER_8_DW_OFFSET    = 12'd384;
localparam [ADDR_PARAM-1:0] LAYER_8_PW_2_OFFSET  = 12'd768;
localparam [ADDR_PARAM-1:0] LAYER_9_PW_1_OFFSET  = 12'd832;
localparam [ADDR_PARAM-1:0] LAYER_9_DW_OFFSET    = 12'd1216;
localparam [ADDR_PARAM-1:0] LAYER_9_PW_2_OFFSET  = 12'd1600;
localparam [ADDR_PARAM-1:0] LAYER_10_PW_1_OFFSET = 12'd1664;
localparam [ADDR_PARAM-1:0] LAYER_10_DW_OFFSET   = 12'd2048;
localparam [ADDR_PARAM-1:0] LAYER_10_PW_2_OFFSET = 12'd2432;
///////////////////////////////////////////////////////
// Registers
///////////////////////////////////////////////////////
reg [14:0]  cnt;
reg         enb_0_q;
reg  [4:0]  dw_cnt;
reg         mean_run, std_run, weight_run, bias_run;
reg  [5:0]  mean_phase, std_phase, weight_phase, bias_phase;


/////////////////////////////////////////////////////// 
// Global counter run window
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        cnt <= 15'd0;
    end else begin
        if (state == PW_1 || state == DW || state == PW_2)
            cnt <= cnt + 15'd1;
        else
            cnt <= 15'd0;
    end
end



///////////////////////////////////////////////////////
// Sequential: address generation
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        enb_0_q        <= 1'b0;
        dw_cnt            <= 5'd0;

        addra_0        <= 0;
        addrb_0        <= 0;
        addra_1        <= 0;
        addrb_1        <= 0;
        addra_w0       <= 0;
        addra_w1       <= 0;
        addra_w2       <= 0;
                          
        addra_biassubam_0   <= 0;
        addra_wdivstd_0     <= 0;

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
                addra_0         <= 0;
                addrb_0         <= 0;
                addra_1         <= 0;
                addrb_1         <= 0;
                addra_w0        <= 0;
                addra_w1        <= 0;
                addra_w2        <= 0;
                if(layer_state == LAYER_8) begin
                    addra_biassubam_0 <= 0;
                    addra_wdivstd_0   <= 0;
                end
                else if(layer_state == LAYER_9) begin
                    addra_biassubam_0 <= LAYER_9_PW_1_OFFSET;
                    addra_wdivstd_0   <= LAYER_9_PW_1_OFFSET;      
                end
                else if(layer_state == LAYER_10) begin
                    addra_biassubam_0 <= LAYER_10_PW_1_OFFSET;
                    addra_wdivstd_0   <= LAYER_10_PW_1_OFFSET;      
                end
            end

            ///////////////////////////////////////////////////////
            // PW_1
            ///////////////////////////////////////////////////////
            PW_1: begin
                // BRAM_A, BRAM_W_0
                if (enb_0) begin
                    addrb_0  <= (addrb_0 >= 6'd63) ? 0 : addrb_0 + 1'b1;
                    if(addra_w0 == 24575) begin
                        addra_w0 <= 0;
                    end
                    else begin
                        addra_w0 <= addra_w0 + 1'b1;
                    end
                end
                else begin
                    addrb_0  <= 0;
                    addra_w0 <= {ADDR_WMEM{1'b0}};
                end

                if(addra_w0 == 15'd24576) begin
                    pw_1_read_done <= 1;
                end
                
                if (save_valid) begin
                    addra_1 <= addra_1 + 1'b1;
                    addra_biassubam_0   <= addra_biassubam_0 + 1;
                    addra_wdivstd_0     <= addra_wdivstd_0 + 1;
                end
            end // PW_1

            ///////////////////////////////////////////////////////
            // PW_1_RST
            ///////////////////////////////////////////////////////
            PW_1_RST: begin
                addra_1             <= 0;      
                
                if(layer_state == LAYER_8) begin
                    addra_biassubam_0 <= LAYER_8_DW_OFFSET;
                    addra_wdivstd_0   <= LAYER_8_DW_OFFSET;
                end
                else if(layer_state == LAYER_9) begin
                    addra_biassubam_0 <= LAYER_9_DW_OFFSET;
                    addra_wdivstd_0   <= LAYER_9_DW_OFFSET;      
                end
                else if(layer_state == LAYER_10) begin
                    addra_biassubam_0 <= LAYER_10_DW_OFFSET;
                    addra_wdivstd_0   <= LAYER_10_DW_OFFSET;      
                end
                
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
                            if(addra_w1 == 3455) begin
                                addra_w1 <= 0;
                            end
                            else begin
                                addra_w1 <= addra_w1 + 1'b1;
                            end
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
                    
                // std: base=46, every 15
                if (cnt == 15'd51) begin
                    addra_wdivstd_0 <= addra_wdivstd_0 + 1'b1; std_run<=1'b1; std_phase<=4'd0;
                end else if (std_run) begin
                    if (std_phase == 6'd28) begin
                        addra_wdivstd_0 <= addra_wdivstd_0 + 1'b1; std_phase <= 4'd0;
                    end else begin
                        std_phase <= std_phase + 1'b1;
                    end
                end
                
                // mean: base=35+1, every 15
                if (cnt == 15'd59) begin
                    addra_biassubam_0 <= addra_biassubam_0 + 1'b1; mean_run<=1'b1; mean_phase<=6'd0;
                end else if (mean_run) begin
                    if (mean_phase == 6'd28) begin
                        addra_biassubam_0 <= addra_biassubam_0 + 1'b1; mean_phase <= 6'd0;
                    end else begin
                        mean_phase <= mean_phase + 6'b1;
                    end
                end


            end // DW

            ///////////////////////////////////////////////////////
            // DW_RST
            ///////////////////////////////////////////////////////
            DW_RST: begin
                addrb_0        <= {INPUT_CHANNEL{1'b0}};
                addra_1        <= {ADDR_CHANNEL{1'b0}};
                
                if(layer_state == LAYER_8) begin
                    addra_biassubam_0 <= LAYER_8_PW_2_OFFSET;
                    addra_wdivstd_0   <= LAYER_8_PW_2_OFFSET;
                end
                else if(layer_state == LAYER_9) begin
                    addra_biassubam_0 <= LAYER_9_PW_2_OFFSET;
                    addra_wdivstd_0   <= LAYER_9_PW_2_OFFSET;      
                end
                else if(layer_state == LAYER_10) begin
                    addra_biassubam_0 <= LAYER_10_PW_2_OFFSET;
                    addra_wdivstd_0   <= LAYER_10_PW_2_OFFSET;      
                end

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
                    if(addra_w2 == 24575) begin 
                        addra_w2 <= 0;
                    end
                    else begin
                        addra_w2 <= addra_w2 + 1'b1;
                    end
                end 
                else begin
                    addrb_1  <= {ADDR_CHANNEL{1'b0}};
                    addra_w2 <= {ADDR_WMEM{1'b0}};
                end
                    
                if(addra_w2 == 15'd24575) begin
                    pw_2_read_done <= 1;
                end
               

                                                                
                // PW2 uses skip_valid
                if (result_save_valid) begin
                    addra_0        <= addra_0 + 1'b1;
                    addra_biassubam_0   <= addra_biassubam_0 + 1;
                    addra_wdivstd_0     <= addra_wdivstd_0 + 1;
                end

                // addrb_0++ on enb_0 falling edge (64 times total)
                if (enb_0_q && !enb_0)
                    addrb_0 <= (addrb_0 == 6'd63) ? 0 : addrb_0 + 1;
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
                dw_read_done   <= 0;
                pw_2_read_done <= 0;
                dw_cnt         <= 0;
                addra_0        <= 0;
                addrb_0        <= 0;
                addra_1        <= 0;
                addrb_1        <= 0;
                addra_w0       <= 0;
                addra_w1       <= 0;
                addra_w2       <= 0;
                addra_biassubam_0   <= 0;
                addra_wdivstd_0     <= 0;
            end
        endcase
    end
end

endmodule