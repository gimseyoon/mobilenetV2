`timescale 1ns / 1ps

module enable_counter #(
    parameter IO_WIDTH      = 18,
    parameter ROW           = 14,
    parameter COLUMN        = 14,
    parameter PIXEL         = ROW * COLUMN,              // 196
    parameter W_WIDTH       = 17,
    parameter ADDR_CHANNEL  = $clog2(384),               // 9
    parameter ADDR_WMEM     = $clog2(384 * 64),          // 15
    parameter ADDR_W1_MEM   = $clog2(384 * 9)            // 12
)(
    input                       clk,
    input                       rst_n,
    input         [2:0]         state,
    input                       save_valid,
    input                       result_save_valid,
    input                       pw_1_read_done,
    input                       dw_read_done,
    input                       pw_2_read_done,
    // BRAM 0
    output                      ena_0,
    output                      wea_0,
    output reg                  enb_0,

    // BRAM 1
    output                      ena_1,
    output                      wea_1,
    output reg                  enb_1,

    // Weight BRAM
    output reg                  ena_w0,
    output reg                  ena_w1,
    output reg                  ena_w2,

    // BN Parameters BRAM
    output reg                  ena_bias_0,
    output reg                  ena_mean_0,
    output reg                  ena_std_0,
    output reg                  ena_weight_0
);

///////////////////////////////////////////////////////
// States
///////////////////////////////////////////////////////
localparam IDLE     = 3'b000,
           PW_1     = 3'b001,
           PW_1_RST = 3'b010,
           DW       = 3'b011,
           DW_RST   = 3'b100,
           PW_2     = 3'b101,
           PW_2_RST = 3'b110,
           EXPORT   = 3'b111;

///////////////////////////////////////////////////////
// Pulse generator params for enb_0 during PW_2
///////////////////////////////////////////////////////
localparam [14:0] START   = 15'd487;   // first trigger
localparam [8:0]  PERIOD  = 9'd384;    // interval
localparam [8:0]  HIGHLEN = 9'd2;      // high length (2 clocks)
localparam [6:0]  REPEAT  = 7'd64;     // total pulses

///////////////////////////////////////////////////////
// Ena/Wea of BRAM ports (combinational)
///////////////////////////////////////////////////////
assign ena_1 = ((state == PW_1) || (state == DW)) && save_valid;
assign wea_1 = ((state == PW_1) || (state == DW)) && save_valid;
assign ena_0 = (state == PW_2) && result_save_valid;
assign wea_0 = (state == PW_2) && result_save_valid;

///////////////////////////////////////////////////////
// Registers for pulse generator
///////////////////////////////////////////////////////

reg [14:0]  cnt;
reg         active;     // PW_2 pulse block active
reg [8:0]   phase;      // 0..PERIOD-1
reg [5:0]   rep_cnt;    // 0..63


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
// Sequential logic
///////////////////////////////////////////////////////

always@(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        ena_bias_0   <= 1'b0;
        ena_mean_0   <= 1'b0;
        ena_std_0    <= 1'b0;
        ena_weight_0 <= 1'b0;
    end
    else begin
        // BN param enable window
        if (state == IDLE) begin
            ena_bias_0   <= 1'b0;
            ena_mean_0   <= 1'b0;
            ena_std_0    <= 1'b0;
            ena_weight_0 <= 1'b0;
        end else begin
            ena_bias_0   <= 1'b1;
            ena_mean_0   <= 1'b1;
            ena_std_0    <= 1'b1;
            ena_weight_0 <= 1'b1;
        end
    end
end

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        enb_0        <= 1'b0;
        enb_1        <= 1'b0;
        ena_w0       <= 1'b0;
        ena_w1       <= 1'b0;
        ena_w2       <= 1'b0;

        active       <= 1'b0;
        phase        <= 9'd0;
        rep_cnt      <= 6'd0;
    end else begin
        case (state)
            ///////////////////////////////////////////////////////
            // IDLE
            ///////////////////////////////////////////////////////
            IDLE: begin
                enb_0        <= 1'b0;
                enb_1        <= 1'b0;
                ena_w0       <= 1'b0;
                ena_w1       <= 1'b0;
                ena_w2       <= 1'b0;

                active       <= 1'b0;
                phase        <= 9'd0;
                rep_cnt      <= 6'd0;
            end

            ///////////////////////////////////////////////////////
            // PW_1
            ///////////////////////////////////////////////////////
            PW_1: begin
                // BRAM_A read + W0 enable
                if (pw_1_read_done) begin
                    enb_0  <= 0;
                    ena_w0 <= 0;
                end else begin
                    enb_0  <= 1;
                    ena_w0 <= 1;
                end
            end

            ///////////////////////////////////////////////////////
            // DW
            ///////////////////////////////////////////////////////
            DW: begin
                // BRAM_B read + W1 enable
                if (dw_read_done) begin
                    enb_1  <= 0;
                    ena_w1 <= 0;
                end
                else begin
                    enb_1  <= 1;
                    ena_w1 <= 1;
                end
            end

            ///////////////////////////////////////////////////////
            // PW_2
            ///////////////////////////////////////////////////////
            PW_2: begin
                // BRAM_B read + W2 enable
                if (pw_2_read_done) begin
                    enb_1  <= 0;
                    ena_w2 <= 0;
                end else begin
                    enb_1  <= 1;
                    ena_w2 <= 1;
                end

                // Two-clock pulses on enb_0, every 384 cycles, 64 times, starting at 473
                if (!active) begin
                    enb_0 <= 1'b0;
                    if (cnt == START) begin
                        active  <= 1'b1;
                        phase   <= 9'd0;
                        rep_cnt <= 6'd0;
                    end
                end 
                else begin
                    enb_0 <= (phase < HIGHLEN);

                    if (phase == (PERIOD - 1)) begin
                        phase <= 9'd0;
                        if (rep_cnt == (REPEAT - 1)) begin
                            active <= 1'b0;  // done (64 pulses)
                            enb_0  <= 1'b0;
                        end else begin
                            rep_cnt <= rep_cnt + 1'b1;
                        end
                    end else begin
                        phase <= phase + 1'b1;
                    end
                end
            
            end
            


            ///////////////////////////////////////////////////////
            // default
            ///////////////////////////////////////////////////////
            default: begin
                enb_0        <= 1'b0;
                enb_1        <= 1'b0;
                ena_w0       <= 1'b0;
                ena_w1       <= 1'b0;
                ena_w2       <= 1'b0;

                active       <= 1'b0;
                phase        <= 9'd0;
                rep_cnt      <= 6'd0;
            end
        endcase
    end
end

endmodule
