`timescale 1ns / 1ps

module accumulator #(
    parameter IO_WIDTH       = 18,
    parameter ROW            = 14,
    parameter COLUMN         = 14,
    parameter PIXEL          = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH        = 17,
    parameter ADDR_CHANNEL   = $clog2(384),               // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM      = $clog2(384 * 64),          // 15 (for 64*384 = 24576)
    parameter ADDR_W1_MEM    = $clog2(384 * 9)            // 12 (for 9*384 = 3456)
)(
    input                                       clk,
    input                                       rst_n,
    input               [2:0]                   state,
    input   signed      [IO_WIDTH*PIXEL-1:0]    mul_out,       // [3528-1:0]
    output  reg                                 bn_en,
    output  reg                                 pw_1_valid,
    output  reg                                 pw_1_done,
    output  reg                                 dw_valid,
    output  reg                                 dw_done,
    output  reg                                 pw_2_valid,
    output  reg                                 pw_2_done,
    output  signed     [IO_WIDTH*PIXEL-1:0]     acc_out        // [3528-1:0]
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
// Registers / wires
///////////////////////////////////////////////////////
reg  [2:0]                         local_state;
reg  [ADDR_CHANNEL-1:0]            channel_num;
reg  [4:0]                         bn_en_cnt;
reg  [ADDR_CHANNEL-1:0]            acc_cnt;
reg  signed [IO_WIDTH-1:0]         acc_out_reg [0:PIXEL-1];
reg  [11:0]                        state_delay;
wire                               pw_1_en;
wire                               dw_en;
wire                               pw_2_en;
reg  signed [IO_WIDTH-1:0]         mul_out_reg [0:PIXEL-1];



/////////////////////////////////////////////////////// 
// local_state
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        local_state <= 0;
    end else begin
        local_state <= state;
    end
end

/////////////////////////////////////////////////////// 
// Enables from state_delay
///////////////////////////////////////////////////////
assign pw_1_en = (local_state == PW_1) ? state_delay[11] : 1'b0;
assign dw_en   = (local_state == DW  ) ? state_delay[9] : 1'b0;
assign pw_2_en = (local_state == PW_2) ? state_delay[10] : 1'b0;

/////////////////////////////////////////////////////// 
// state_delay (BRAM read wait)
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state_delay <= 0;
    end else begin
        state_delay[11]<= state_delay[10];
        state_delay[10]<= state_delay[9];
        state_delay[9] <= state_delay[8];
        state_delay[8] <= state_delay[7];
        state_delay[7] <= state_delay[6];
        state_delay[6] <= state_delay[5];
        state_delay[5] <= state_delay[4];
        state_delay[4] <= state_delay[3];
        state_delay[3] <= state_delay[2];
        state_delay[2] <= state_delay[1];
        state_delay[1] <= state_delay[0];
        if (local_state == PW_1 || local_state == DW || local_state == PW_2)
            state_delay[0] <= 1;
        else
            state_delay     <= 0;
    end
end

/////////////////////////////////////////////////////// 
// Unpack mul_out -> mul_out_reg[PIXEL]
///////////////////////////////////////////////////////
integer p;

always @(posedge clk or negedge rst_n) begin
  if (!rst_n) begin
    for (p = 0; p < PIXEL; p = p + 1)
      mul_out_reg[p] <= 0;
  end else begin
    for (p = 0; p < PIXEL; p = p + 1)
      mul_out_reg[p] <= mul_out[p*IO_WIDTH +: IO_WIDTH];
  end
end
/////////////////////////////////////////////////////// 
// Pack acc_out_reg[] -> acc_out
///////////////////////////////////////////////////////
genvar i;
generate
    for (i = 0; i < PIXEL; i = i + 1) begin : PACK_OUTPUT
        assign acc_out[IO_WIDTH*(i+1)-1 : IO_WIDTH*i] = acc_out_reg[i];
    end
endgenerate

/////////////////////////////////////////////////////// 
// Main accumulate / valid / done / bn_en
///////////////////////////////////////////////////////
integer k, x, y;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        for (k = 0; k < PIXEL; k = k + 1) acc_out_reg[k] <= 0;
        bn_en_cnt   <= 0;
        acc_cnt     <= 0;
        channel_num <= 0;
        bn_en       <= 0;
        pw_1_valid  <= 0; pw_1_done <= 0;
        dw_valid    <= 0; dw_done   <= 0;
        pw_2_valid  <= 0; pw_2_done <= 0;
    end else begin
        case (local_state)
            IDLE: begin
                for (k = 0; k < PIXEL; k = k + 1) acc_out_reg[k] <= 0;
                bn_en_cnt   <= 0;
                acc_cnt     <= 0;
                channel_num <= 0;
                bn_en       <= 0;
                pw_1_valid  <= 0; pw_1_done <= 0;
                dw_valid    <= 0; dw_done   <= 0;
                pw_2_valid  <= 0; pw_2_done <= 0;
            end

            PW_1: begin

                for (k = 0; k < PIXEL; k = k + 1) begin
                    if (acc_cnt >= 9'd63)  acc_out_reg[k] <= mul_out_reg[k];      
                    else                   acc_out_reg[k] <= acc_out_reg[k] + mul_out_reg[k];
                end


                if (pw_1_en) begin
                    acc_cnt     <= (acc_cnt == 9'd63) ? 0 : acc_cnt + 1'b1;
                    pw_1_valid  <= (channel_num != 9'd511) ? (acc_cnt == 9'd62) : 0;
                    pw_1_done   <= (acc_cnt == 9'd62) && (channel_num == 9'd383);

                    if (pw_1_valid) begin
                        if       (channel_num == 9'd383) channel_num <= 11'd511;
                        else if (channel_num == 9'd511) channel_num <= channel_num;
                        else                             channel_num <= channel_num + 1'b1;
                    end
                end else begin
                    acc_cnt     <= 0;
                    pw_1_valid  <= 0;
                    pw_1_done   <= 0;
                    channel_num <= 0;
                end

                if (pw_1_valid) begin
                    bn_en <= (channel_num != 9'd511);
                end else if (bn_en_cnt == 5'd27) begin
                    bn_en <= 1'b0;
                end

                if (bn_en) begin
                    if (bn_en_cnt == 5'd27) bn_en_cnt <= 0;
                    else                    bn_en_cnt <= bn_en_cnt + 1'b1;
                end else begin
                    bn_en_cnt <= 0;
                end
            end // PW_1

            DW: begin
                    case (acc_cnt)
                        // 1) DR(13??13) -> (y+1,x+1) overwrite
                        4'd1: begin
                            for (y = 0; y < 13; y = y + 1)
                                for (x = 0; x < 13; x = x + 1)
                                    acc_out_reg[(y+1)*14 + (x+1)] <= mul_out_reg[y*14 + x];
                        end
                        // 2) D(14??13) -> (y+1,x) +=
                        4'd2: begin
                            for (y = 0; y < 13; y = y + 1)
                                for (x = 0; x < 14; x = x + 1)
                                    acc_out_reg[(y+1)*14 + x] <= acc_out_reg[(y+1)*14 + x] + mul_out_reg[y*14 + x];
                        end
                        // 3) DL(13??13) -> (y+1,x-1) +=
                        4'd3: begin
                            for (y = 0; y < 13; y = y + 1)
                                for (x = 1; x < 14; x = x + 1)
                                    acc_out_reg[(y+1)*14 + (x-1)] <= acc_out_reg[(y+1)*14 + (x-1)] + mul_out_reg[y*14 + x];
                        end
                        // 4) R(14??13) -> (y,x+1) +=
                        4'd4: begin
                            for (y = 0; y < 14; y = y + 1)
                                for (x = 0; x < 13; x = x + 1)
                                    acc_out_reg[y*14 + (x+1)] <= acc_out_reg[y*14 + (x+1)] + mul_out_reg[y*14 + x];
                        end
                        // 5) ALL(14??14) +=
                        4'd5: begin
                            for (p = 0; p < PIXEL; p = p + 1)
                                acc_out_reg[p] <= acc_out_reg[p] + mul_out_reg[p];
                        end
                        // 6) L(14??13) -> (y,x-1) +=
                        4'd6: begin
                            for (y = 0; y < 14; y = y + 1)
                                for (x = 1; x < 14; x = x + 1)
                                    acc_out_reg[y*14 + (x-1)] <= acc_out_reg[y*14 + (x-1)] + mul_out_reg[y*14 + x];
                        end
                        // 7) UR(13??13) -> (y-1,x+1) +=
                        4'd7: begin
                            for (y = 1; y < 14; y = y + 1)
                                for (x = 0; x < 13; x = x + 1)
                                    acc_out_reg[(y-1)*14 + (x+1)] <= acc_out_reg[(y-1)*14 + (x+1)] + mul_out_reg[y*14 + x];
                        end
                        // 8) U(13??14) -> (y-1,x) +=
                        4'd8: begin
                            for (y = 1; y < 14; y = y + 1)
                                for (x = 0; x < 14; x = x + 1)
                                    acc_out_reg[(y-1)*14 + x] <= acc_out_reg[(y-1)*14 + x] + mul_out_reg[y*14 + x];
                        end
                        // 9) UL(13??13) -> (y-1,x-1) +=
                        4'd9: begin
                            for (y = 1; y < 14; y = y + 1)
                                for (x = 1; x < 14; x = x + 1)
                                    acc_out_reg[(y-1)*14 + (x-1)] <= acc_out_reg[(y-1)*14 + (x-1)] + mul_out_reg[y*14 + x];
                        end
                        default: begin
                            for (k = 0; k < PIXEL; k = k + 1) acc_out_reg[k] <= 0;
                        end
                    endcase

                    if (channel_num == 9'd511 && acc_cnt == 9'd9) begin
                        bn_en <= 1'b0;
                    end else if (channel_num != 9'd511 && dw_valid) begin
                        bn_en <= 1'b1;
                    end else if (bn_en_cnt == 5'd27) begin
                        bn_en <= 1'b0;
                    end

                    if (bn_en) begin
                        if (bn_en_cnt == 5'd27) bn_en_cnt <= 0;
                        else                    bn_en_cnt <= bn_en_cnt + 1'b1;
                    end else begin
                        bn_en_cnt <= 0;
                    end

                    if (dw_en) begin
                        acc_cnt   <= (acc_cnt == 9'd29) ? 5'd1 : acc_cnt + 1'b1;
                        dw_valid  <= (acc_cnt == 9'd9);
                        dw_done   <= (acc_cnt == 9'd9) && (channel_num == 9'd383);

                        if (acc_cnt == 9'd29) begin
                            if      (channel_num == 9'd383) channel_num <= 9'd511;
                            else if (channel_num == 9'd511) channel_num <= channel_num;
                            else                             channel_num <= channel_num + 1'b1;
                        end
                    end else begin
                        acc_cnt   <= 0;
                        dw_valid  <= 0;
                        dw_done   <= 0;
                        channel_num <= 0;
                    end
            end // DW
            
            
            PW_2: begin

                for (k = 0; k < PIXEL; k = k + 1) begin
                    if (acc_cnt >= 9'd383)  acc_out_reg[k] <= mul_out_reg[k];            
                    else                    acc_out_reg[k] <= acc_out_reg[k] + mul_out_reg[k];
                end

                if (pw_2_en) begin
                    acc_cnt    <= (acc_cnt == 9'd383) ? 0 : acc_cnt + 1'b1;
                    pw_2_valid <= (acc_cnt == 9'd382);
                    pw_2_done  <= (acc_cnt == 9'd382) && (channel_num == 9'd63);

                    if (pw_2_valid) begin
                        if      (channel_num == 9'd63)  channel_num <= 9'd511;
                        else if (channel_num == 9'd511) channel_num <= channel_num;
                        else                              channel_num <= channel_num + 1'b1;
                    end
                end else begin
                    acc_cnt    <= 511;
                    pw_2_valid <= 0;
                    pw_2_done  <= 0;
                    channel_num <= 0;
                end

                if (pw_2_valid) begin
                    bn_en <= (channel_num != 9'd511);
                end else if (bn_en_cnt == 5'd27) begin
                    bn_en <= 1'b0;
                end

                if (bn_en) begin
                    if (bn_en_cnt == 5'd27) bn_en_cnt <= 0;
                    else                    bn_en_cnt <= bn_en_cnt + 1'b1;
                end else begin
                    bn_en_cnt <= 0;
                end
            end

            EXPORT: begin
                // reserved
            end

            default: begin
                bn_en_cnt   <= 0;
                acc_cnt     <= 0;
                channel_num <= 9'd511;
                bn_en       <= 0;
                pw_1_valid  <= 0; pw_1_done <= 0;
                dw_valid    <= 0; dw_done   <= 0;
                pw_2_valid  <= 0; pw_2_done <= 0;
            end
        endcase
    end
end

endmodule