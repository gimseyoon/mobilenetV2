`timescale 1ns / 1ps

module mobilenetV2 #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),         // 8 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64)       // 15 (for 64*384 = 24576)
)(
    input clk,
    input rst_n,
    input start
    );
    
    
    
//////////////////////////////////////////////////
    wire [2:0] state;
    wire bram_select;
    
//////////////////////////////////////////////////
// multiplier
    wire signed [IO_WIDTH * PIXEL - 1 : 0] mul_in;   // [3920-1 : 0], 3920 bit
    wire signed [W_WIDTH - 1 : 0] mul_weight;        // [17-1:0], 17 bit
    wire signed [IO_WIDTH * PIXEL - 1 : 0] mul_out; // [3920-1 : 0], 3920 bit
    
//////////////////////////////////////////////////
//accumultaor
    wire signed [IO_WIDTH * PIXEL - 1 : 0] acc_out;
    wire [ADDR_CHANNEL-1 : 0] channel_num;
    wire pw_1_done;
    wire bn_relu_valid;
    wire [ADDR_CHANNEL -1 :0] acc_cnt;
//////////////////////////////////////////////////
// BN_RELU
    wire signed [IO_WIDTH * PIXEL - 1 : 0] bn_relu_out;
    wire save_valid;
    
//////////////////////////////////////////////////
// bram_0
    wire                                        ena_0;
    wire                                        wea_0;
    wire  [ADDR_CHANNEL-1:0]                    addra_0;
    wire  signed [IO_WIDTH * PIXEL - 1 : 0]     dina_0;
    wire                                        enb_0;
    wire  [ADDR_CHANNEL-1:0]                    addrb_0;
    wire  signed [IO_WIDTH * PIXEL - 1 : 0]     doutb_0;
// bram_1
    wire                                        ena_1;
    wire                                        wea_1;
    wire  [ADDR_CHANNEL-1:0]                    addra_1;
    wire  signed [IO_WIDTH * PIXEL - 1 : 0]     dina_1;
    wire                                        enb_1;
    wire  [ADDR_CHANNEL-1:0]                    addrb_1;
    wire  signed [IO_WIDTH * PIXEL - 1 : 0]     doutb_1;
// bram_w
    wire                                        ena_w0;
    wire  [ADDR_WMEM-1:0]                       addra_w0; 
    wire  signed [W_WIDTH-1:0]                  douta_w0;
// bram_bias
    wire                                        ena_bias_0;
    wire  [ADDR_CHANNEL-1 : 0]                  addra_bias_0;
    wire  [31:0]                                douta_bias_0;
// bram_mean
    wire                                        ena_mean_0;
    wire  [ADDR_CHANNEL-1 : 0]                  addra_mean_0;
    wire  [31:0]                                douta_mean_0;
// bram_std
    wire                                        ena_std_0;
    wire  [ADDR_CHANNEL-1 : 0]                  addra_std_0;
    wire  [31:0]                                douta_std_0;
// bram_weight
    wire                                        ena_weight_0;
    wire  [ADDR_CHANNEL-1 : 0]                  addra_weight_0;
    wire  [31:0]                                douta_weight_0;
//////////////////////////////////////////////////

    assign mul_weight = douta_w0;
    assign mul_in = doutb_0;
    


//////////////////////////////////////////////////
// Instantiate FSM
FSM FSM_0 (
    .clk                (clk),
    .rst_n              (rst_n),
    .start              (start),
    .pw_1_done          (pw_1_done),
    .pw_1_bn_relu_done  (pw_1_bn_relu_done),
    .dw_done            (dw_done),
    .dw_bn_relu_done    (dw_bn_relu_done),
    .pw_2_done          (pw_2_done),
    .pw_2_bn_done       (pw_2_bn_done),
    .layer_done         (layer_done),  
    .bram_select        (bram_select),
    .state              (state)
);
    

//////////////////////////////////////////////////
// Instantiate addr_counter

glbl_ctrl glbl_ctrl_0 (
    .clk                (clk),
    .rst_n              (rst_n),
    .state              (state),
    .save_valid         (save_valid),
    .channel_num        (channel_num),
    .bram_select        (bram_select),
    .acc_out            (acc_out),
    .acc_cnt            (acc_cnt),
// BRAM A
    .ena_0              (ena_0),
    .wea_0              (wea_0),
    .addra_0            (addra_0),
    .dina_0             (dina_0),
    .enb_0              (enb_0),
    .addrb_0            (addrb_0),
// BRAM B
    .ena_1              (ena_1),
    .wea_1              (wea_1),
    .addra_1            (addra_1),
    .dina_1             (dina_1),    
    .enb_1              (enb_1),
    .addrb_1            (addrb_1),
// BRAM W
    .ena_w0             (ena_w0),
    .addra_w0           (addra_w0),
    
// BRAM bias
    .ena_bias_0         (ena_bias_0),
    .addra_bias_0       (addra_bias_0),
// BRAM mean
    .ena_mean_0         (ena_mean_0),
    .addra_mean_0       (addra_mean_0),
// BRAM std
    .ena_std_0          (ena_std_0),
    .addra_std_0        (addra_std_0),
// BRAM weight
    .ena_weight_0       (ena_weight_0),
    .addra_weight_0     (addra_weight_0)
);



    
//////////////////////////////////////////////////
// Instantiate accumulator
accumulator accumulator_0 (
    .clk                (clk),
    .rst_n              (rst_n),
    .state              (state),
    .mul_out            (mul_out),
    .acc_cnt            (acc_cnt),
    .bn_relu_valid      (bn_relu_valid),
    .channel_num        (channel_num),
    .pw_1_done          (pw_1_done),
    .acc_out            (acc_out)
);
    
//////////////////////////////////////////////////
// Instantiate multiplier
multiplier multiplier_0 (
    .clk                (clk),
    .rst_n              (rst_n),
    .mul_in             (mul_in),
    .mul_weight         (mul_weight),
    .mul_out            (mul_out)
);
    
//////////////////////////////////////////////////
// Instantiate BN_RELU
BN_RELU BN_RELU_0 (
    .clk                (clk),
    .rst_n              (rst_n),
    .bn_relu_valid      (bn_relu_valid),
    .mean               (douta_mean_0),
    .weight             (douta_weight_0),
    .bias               (douta_bias_0),
    .std                (douta_std_0),
    .acc_out            (acc_out), 
    .bn_relu_out        (bn_relu_out), 
    .save_valid         (save_valid)
);

    
//////////////////////////////////////////////////
// BRAM (INPUT, OUTPUT, WEIGHT)

blk_mem_gen_0 IO_bram_A (
  .clka                 (clk),    // input wire clka
  .ena                  (ena_0),      // input wire ena
  .wea                  (wea_0),      // input wire [0 : 0] wea
  .addra                (addra_0),  // input wire [8 : 0] addra
  .dina                 (dina_0),    // input wire [3527 : 0] dina
  .clkb                 (clk),    // input wire clkb
  .enb                  (enb_0),      // input wire enb
  .addrb                (addrb_0),  // input wire [8 : 0] addrb
  .doutb                (doutb_0)  // output wire [3527 : 0] doutb
);


blk_mem_gen_1 IO_bram_B (
  .clka                 (clk),    // input wire clka
  .ena                  (ena_1),      // input wire ena
  .wea                  (wea_1),      // input wire [0 : 0] wea
  .addra                (addra_1),  // input wire [8 : 0] addra
  .dina                 (dina_1),    // input wire [3527 : 0] dina
  .clkb                 (clk),    // input wire clkb
  .enb                  (enb_1),      // input wire enb
  .addrb                (addrb_1),  // input wire [8 : 0] addrb
  .doutb                (doutb_1)  // output wire [3527 : 0] doutb
);


blk_mem_gen_2 W_bram_w0 (
  .clka                 (clk),    // input wire clka
  .ena                  (ena_w0),      // input wire ena
  .addra                (addra_w0),  // input wire [14 : 0] addra
  .douta                (douta_w0)  // output wire [16 : 0] douta
);

///////////////////////////////////////////////////////////////
// bram_PARAMETER 
bram_bias bram_bias_0 (
  .clka(clk),    // input wire clka
  .ena(ena_bias_0),      // input wire ena
  .addra(addra_bias_0),  // input wire [8 : 0] addra
  .douta(douta_bias_0)  // output wire [31 : 0] douta
);

bram_mean bram_mean_0 (
  .clka(clk),    // input wire clka
  .ena(ena_mean_0),      // input wire ena
  .addra(addra_mean_0),  // input wire [8 : 0] addra
  .douta(douta_mean_0)  // output wire [31 : 0] douta
);

bram_std bram_std_0 (
  .clka(clk),    // input wire clka
  .ena(ena_std_0),      // input wire ena
  .addra(addra_std_0),  // input wire [8 : 0] addra
  .douta(douta_std_0)  // output wire [31 : 0] douta
);

bram_weight bram_weight_0 (
  .clka(clk),    // input wire clka
  .ena(ena_weight_0),      // input wire ena
  .addra(addra_weight_0),  // input wire [8 : 0] addra
  .douta(douta_weight_0)  // output wire [31 : 0] douta
);

/////////////////////////////////////////////////////////////

endmodule