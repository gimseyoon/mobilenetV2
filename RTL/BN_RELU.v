`timescale 1ns / 1ps

module BN_RELU #(
    parameter IO_WIDTH      = 18,
    parameter ROW           = 14,
    parameter COLUMN        = 14,
    parameter PIXEL         = ROW * COLUMN,                 // 14 * 14 = 196
    parameter W_WIDTH       = 17,
    parameter ADDR_CHANNEL  = $clog2(384),                  // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM     = $clog2(384 * 64)              // 15 (for 64*384 = 24576)
)(
    input                                   clk,
    input                                   rst_n,
    input           [2:0]                   state,
    input                                   pw_1_valid,
    input                                   dw_valid,
    input                                   pw_2_valid,
    input                                   bn_en,
    input   signed  [31:0]                  biassubam,
    input   signed  [31:0]                  wdivstd,
    input   signed  [IO_WIDTH*PIXEL-1:0]    acc_out,        // [3528-1 : 0]
    output  signed  [IO_WIDTH*PIXEL-1:0]    bn_relu_out,    // [3528-1 : 0]
    output  reg                             save_valid,
    output  reg                             skip_valid,
    output  reg                             pw_1_bn_relu_done,
    output  reg                             dw_bn_relu_done,
    output  reg                             pw_2_bn_done
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
// Registers / Wires
///////////////////////////////////////////////////////
reg          [2:0]                  local_state;
reg          [4:0]                  bn_cnt;             // 0~13
reg          [4:0]                  bn_save_cnt;        // 0~27
reg          [ADDR_CHANNEL-1:0]     bn_channel_num;
reg   signed [IO_WIDTH-1:0]         acc_out_reg       [0:PIXEL-1];
wire  signed [IO_WIDTH-1:0]         acc_selected      [0:6];
wire  signed [IO_WIDTH-1:0]         bn_single_out     [0:6];
reg   signed [IO_WIDTH-1:0]         bn_relu_out_array [0:PIXEL-1];
wire         [6:0]                  valid_single;
wire                                bn_valid;


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
// All 14 lanes valid
///////////////////////////////////////////////////////
assign bn_valid = &valid_single;

/////////////////////////////////////////////////////// 
// Counters, save/skip valid
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        bn_cnt          <= 5'd0;
        bn_save_cnt     <= 5'd0;
        save_valid      <= 1'b0;
        skip_valid      <= 1'b0;
    end else begin
        // bn_cnt : bn_en ???? 0..27 ???
        if (bn_cnt == 5'd27)      bn_cnt <= 5'd0;
        else if (bn_en)           bn_cnt <= bn_cnt + 1'b1;

        // bn_save_cnt : bn_valid?? ???? 0..13 ???, ???? 0
        if (bn_valid) begin
            if (bn_save_cnt == 5'd27) bn_save_cnt <= 5'd0;
            else                       bn_save_cnt <= bn_save_cnt + 1'b1;
        end else begin
            bn_save_cnt <= 5'd0;
        end

        // save_valid / skip_valid
        if (local_state == PW_1 || local_state == DW) begin
            if (bn_valid && bn_save_cnt == 5'd27) save_valid <= 1'b1;
            else                                  save_valid <= 1'b0;
        end 
        else begin
            save_valid <= 1'b0;
        end

        if (local_state == PW_2) begin
            if (bn_valid && bn_save_cnt == 5'd27) skip_valid <= 1'b1;
            else                                  skip_valid <= 1'b0;
        end 
        else begin
            skip_valid <= 1'b0;
        end
    end
end

/////////////////////////////////////////////////////// 
// Capture acc_out into registers on valid
///////////////////////////////////////////////////////
genvar i;
generate
    for (i = 0; i < PIXEL; i = i + 1) begin : REG_ARRAY_ASSIGN
        always @(posedge clk or negedge rst_n) begin
            if (!rst_n) begin
                acc_out_reg[i] <= {IO_WIDTH{1'b0}};
            end 
            else  begin
                if (pw_1_valid || dw_valid || pw_2_valid) begin
                    acc_out_reg[i] <= acc_out[IO_WIDTH*(i+1)-1 : IO_WIDTH*i];
                end 
            end
        end //always
    end
endgenerate

/////////////////////////////////////////////////////// 
// Select 7 elements (current row)
///////////////////////////////////////////////////////
genvar s;
generate
    for (s = 0; s < 7; s = s + 1) begin : SELECT_ACC
        assign acc_selected[s] = acc_out_reg[bn_cnt*7 + s];
    end
endgenerate

/////////////////////////////////////////////////////// 
// Output pack
///////////////////////////////////////////////////////
genvar k;
generate
    for (k = 0; k < PIXEL; k = k + 1) begin : PACK_OUTPUT
        assign bn_relu_out[IO_WIDTH*(k+1)-1 : IO_WIDTH*k] = bn_relu_out_array[k];
    end
endgenerate

/////////////////////////////////////////////////////// 
// Channel number (0..383) : ?? ?? ???? ?????? +1
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        bn_channel_num <= {ADDR_CHANNEL{1'b0}};
    end else begin
        if (bn_save_cnt == 5'd27) begin
            if (bn_channel_num == 9'd383) bn_channel_num <= {ADDR_CHANNEL{1'b0}};
            else                          bn_channel_num <= bn_channel_num + 1'b1;
        end
    end
end

/////////////////////////////////////////////////////// 
// Writeback into bn_relu_out_array
///////////////////////////////////////////////////////
integer q, t, idx;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        for (q = 0; q < PIXEL; q = q + 1)
            bn_relu_out_array[q] <= {IO_WIDTH{1'b0}};
    end else begin
        if (bn_valid) begin
            for (t = 0; t < 7; t = t + 1) begin
                idx = bn_save_cnt*7 + t;
                if (local_state == PW_1 || local_state == DW) begin
                    // ReLU
                    if (bn_single_out[t][IO_WIDTH-1]) bn_relu_out_array[idx] <= {IO_WIDTH{1'b0}};
                    else                               bn_relu_out_array[idx] <= bn_single_out[t];
                end else if (local_state == PW_2) begin
                    bn_relu_out_array[idx] <= bn_single_out[t];
                end else begin
                    bn_relu_out_array[idx] <= {IO_WIDTH{1'b0}};
                end
            end
        end
    end
end

/////////////////////////////////////////////////////// 
// Done flags
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        pw_1_bn_relu_done <= 1'b0;
        dw_bn_relu_done   <= 1'b0;
        pw_2_bn_done      <= 1'b0;
    end else begin
        pw_1_bn_relu_done <= (local_state == PW_1) && (bn_save_cnt == 5'd27) && (bn_channel_num == 9'd383);
        dw_bn_relu_done   <= (local_state == DW)   && (bn_save_cnt == 5'd27) && (bn_channel_num == 9'd383);
        pw_2_bn_done      <= (local_state == PW_2) && (bn_save_cnt == 5'd27) && (bn_channel_num == 9'd063);
    end
end

/////////////////////////////////////////////////////// 
// 7-way BN units
///////////////////////////////////////////////////////
genvar p;
generate
    for (p = 0; p < 7; p = p + 1) begin : BN_UNIT
        NEW_BN_RELU_SINGLE #(.IO_WIDTH(IO_WIDTH)) u_bn (
            .clk       (clk),
            .bn_en     (bn_en),
            .biassubam (biassubam),
            .wdivstd   (wdivstd),
            .acc_in    (acc_selected[p]),
            .bn_out    (bn_single_out[p]),
            .valid_out (valid_single[p])
        );
    end
endgenerate



endmodule