`timescale 1ns / 1ps

module glbl_ctrl #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter INPUT_CHANNEL = 64,
    parameter ADDR_PARAM = 12,
    parameter ADDR_IN = $clog2(64),
    parameter ADDR_CHANNEL  = $clog2(384),        // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(24576*3),       // 15 (for 64*384 = 24576)
    parameter ADDR_W1_MEM = $clog2(3456*3)       // 12 (for 9*384 = 3456)
)(
    input                                   clk,
    input                                   rst_n,
    input           [2:0]                   state,
    input           [1:0]                   layer_state,
    input                                   save_valid,
    input                                   skip_valid,
    input   signed  [IO_WIDTH*PIXEL-1:0]    acc_out,
    input                                   pw_1_done,
    input                                   result_save_valid,

    /////////////////////////////////////////////////////// 
    // BRAM_A (input/output feature map)
    ///////////////////////////////////////////////////////
    output                                  ena_0,
    output                                  wea_0,
    output             [INPUT_CHANNEL-1:0]  addra_0,
    output reg signed [IO_WIDTH*PIXEL-1:0]  dina_0,
    output                                  enb_0,
    output             [INPUT_CHANNEL-1:0]  addrb_0,

    /////////////////////////////////////////////////////// 
    // BRAM_B (intermediate feature map)
    ///////////////////////////////////////////////////////
    output                                  ena_1,
    output                                  wea_1,
    output              [ADDR_CHANNEL-1:0]  addra_1,
    output reg signed [IO_WIDTH*PIXEL-1:0]  dina_1,
    output                                  enb_1,
    output              [ADDR_CHANNEL-1:0]  addrb_1,

    /////////////////////////////////////////////////////// 
    // Weights BRAMs
    ///////////////////////////////////////////////////////
    output                                  ena_w0,
    output          [ADDR_WMEM-1:0]         addra_w0,
    output                                  ena_w1,
    output          [ADDR_W1_MEM-1:0]       addra_w1,
    output                                  ena_w2,
    output          [ADDR_WMEM-1:0]         addra_w2,

    /////////////////////////////////////////////////////// 
    // BN parameter BRAMs
    ///////////////////////////////////////////////////////
    output                                  ena_biassubam_0,
    output          [ADDR_PARAM-1:0]        addra_biassubam_0,
    output                                  ena_wdivstd_0,
    output          [ADDR_PARAM-1:0]        addra_wdivstd_0
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
           SK       = 3'b111;
           
localparam READY    = 2'b00,
           LAYER_8  = 2'b01,
           LAYER_9  = 2'b10,
           LAYER_10 = 2'b11;
/////////////////////////////////////////////////////// 
// Global counter (0..32767)
///////////////////////////////////////////////////////
reg  [2:0] local_state;
reg  [1:0] local_layer_state;
wire pw_1_read_done;
wire dw_read_done;
wire pw_2_read_done;


/////////////////////////////////////////////////////// 
// local_state
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        local_state <= 0;
        local_layer_state <= 0;
    end else begin
        local_state <= state;
        local_layer_state <= layer_state;
    end
end




/////////////////////////////////////////////////////// 
// addr_counter
///////////////////////////////////////////////////////
addr_counter addr_counter_0 (
    .clk                (clk),
    .rst_n              (rst_n),
    .state              (local_state),
    .layer_state        (local_layer_state),
    .save_valid         (save_valid),
    .skip_valid         (skip_valid),
    .result_save_valid (result_save_valid),
    .enb_0              (enb_0),
    .enb_1              (enb_1),

    .pw_1_read_done     (pw_1_read_done),
    .dw_read_done       (dw_read_done),
    .pw_2_read_done     (pw_2_read_done),
    
    .addra_0            (addra_0),
    .addrb_0            (addrb_0),
    .addra_1            (addra_1),
    .addrb_1            (addrb_1),
    .addra_w0           (addra_w0),
    .addra_w1           (addra_w1),
    .addra_w2           (addra_w2),
    .addra_biassubam_0  (addra_biassubam_0),
    .addra_wdivstd_0    (addra_wdivstd_0)
);
/////////////////////////////////////////////////////// 
// enable_counter
///////////////////////////////////////////////////////
enable_counter enable_counter_0 (
    .clk                (clk),
    .rst_n              (rst_n),
    .state              (local_state),
    .save_valid         (save_valid),
    .result_save_valid  (result_save_valid),
    .pw_1_read_done     (pw_1_read_done),
    .dw_read_done       (dw_read_done),
    .pw_2_read_done     (pw_2_read_done),

    .ena_0              (ena_0),
    .wea_0              (wea_0),
    .enb_0              (enb_0),

    .ena_1              (ena_1),
    .wea_1              (wea_1),
    .enb_1              (enb_1),

    .ena_w0             (ena_w0),
    .ena_w1             (ena_w1),
    .ena_w2             (ena_w2),

    .ena_biassubam_0  (ena_biassubam_0),
    .ena_wdivstd_0    (ena_wdivstd_0)
);

endmodule