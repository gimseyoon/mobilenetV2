`timescale 1ns / 1ps

module mobilenetV2 #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),        // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64),       // 15 (for 64*384 = 24576)
    parameter ADDR_W1_MEM = $clog2(384 * 9)       // 12 (for 9*384 = 3456)
)(
    input clk,
    input rst_n,
    input start
);
    
//////////////////////////////////////////////////

    localparam IDLE         = 3'b000,
               PW_1         = 3'b001,
               PW_1_BN_RELU = 3'b010,
               DW           = 3'b011,
               DW_BN_RELU   = 3'b100,
               PW_2         = 3'b101,
               PW_2_BN      = 3'b110,
               SK           = 3'b111;

//////////////////////////////////////////////////
    
//////////////////////////////////////////////////
// FSM
    wire [2:0] state;
    wire pw_1_bn_relu_done;
    wire dw_bn_relu_done;

    wire pw_2_bn_done;
    wire layer_done;
    wire bn_relu_valid;
//////////////////////////////////////////////////
// multiplier
    reg signed [IO_WIDTH * PIXEL - 1 : 0] mul_in;   // [3920-1 : 0], 3920 bit
    reg signed [W_WIDTH - 1 : 0] mul_weight;        // [17-1:0], 17 bit
    wire signed [IO_WIDTH * PIXEL - 1 : 0] mul_out; // [3920-1 : 0], 3920 bit
    
//////////////////////////////////////////////////
//accumultaor
    wire signed [IO_WIDTH * PIXEL - 1 : 0] acc_out;
    wire [ADDR_CHANNEL-1 : 0] channel_num;
    wire pw_1_valid;   
    wire pw_1_done;
    wire dw_valid;
    wire dw_done;
    wire pw_2_valid;
    wire pw_2_done;
    wire [ADDR_CHANNEL -1 :0] acc_cnt;
    wire save_valid;
//////////////////////////////////////////////////
// BN_RELU
    wire signed [IO_WIDTH * PIXEL - 1 : 0] bn_relu_out;
    wire [3:0] bn_cnt;
    
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
    
    
// bram_W0
    wire                                        ena_w0;
    wire  [ADDR_WMEM-1:0]                       addra_w0; 
    wire  signed [W_WIDTH-1:0]                  douta_w0;
// bram_W1
    wire                                        ena_w1;
    wire  [ADDR_W1_MEM-1:0]                     addra_w1; 
    wire  signed [W_WIDTH-1:0]                  douta_w1;
// bram_W2
    wire                                        ena_w2;
    wire  [ADDR_WMEM-1:0]                       addra_w2; 
    wire  signed [W_WIDTH-1:0]                  douta_w2;
    
    
// bram_bias
    wire                                        ena_bias_0;
    wire  [ADDR_CHANNEL-1 : 0]                  addra_bias_0;
    wire  signed        [31:0]                  douta_bias_0;
// bram_mean
    wire                                        ena_mean_0;
    wire  [ADDR_CHANNEL-1 : 0]                  addra_mean_0;
    wire  signed        [31:0]                  douta_mean_0;
// bram_std
    wire                                        ena_std_0;
    wire  [ADDR_CHANNEL-1 : 0]                  addra_std_0;
    wire  signed        [31:0]                  douta_std_0;
// bram_weight
    wire                                        ena_weight_0;
    wire  [ADDR_CHANNEL-1 : 0]                  addra_weight_0;
    wire  signed        [31:0]                  douta_weight_0;
    
//////////////////////////////////////////////////

    always@(*) begin
        case(state)
            IDLE: begin
                mul_in = 0;
                mul_weight = 0;
            end
            PW_1: begin
                mul_in = doutb_0;
                mul_weight = douta_w0;
            end
            PW_1_BN_RELU: begin
                mul_in = 0;
                mul_weight = 0;
            end
            DW: begin
                mul_in = doutb_1;
                mul_weight = douta_w1;
            end
            default: begin
                mul_in = 0;
                mul_weight = 0;
            end
        endcase
    end    
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
    .state              (state)
);
    

//////////////////////////////////////////////////
// Instantiate addr_counter

glbl_ctrl glbl_ctrl_0 (
    .clk                (clk),
    .rst_n              (rst_n),
    .state              (state),
    .bn_cnt             (bn_cnt),
    .pw_1_valid         (pw_1_valid),   
    .dw_valid           (dw_valid), 
    .save_valid         (save_valid),
    .channel_num        (channel_num),
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
    
    
// BRAM W0
    .ena_w0             (ena_w0),
    .addra_w0           (addra_w0),
// BRAM W1
    .ena_w1             (ena_w1),
    .addra_w1           (addra_w1),
// BRAM W2
    .ena_w2             (ena_w2),
    .addra_w2           (addra_w2),
    
    
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
    .bn_en              (bn_en),
    .acc_cnt            (acc_cnt),
    .pw_1_valid         (pw_1_valid),
    .pw_1_done          (pw_1_done),
    .dw_valid           (dw_valid),
    .dw_done            (dw_done),
    .pw_2_valid         (pw_2_valid),
    .pw_2_done          (pw_2_done),
    .channel_num        (channel_num),
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
    .pw_1_valid         (pw_1_valid),
    .bn_en              (bn_en),
    .mean               (douta_mean_0),
    .weight             (douta_weight_0),
    .bias               (douta_bias_0),
    .std                (douta_std_0),
    .acc_out            (acc_out), 
    .bn_cnt             (bn_cnt),
    .bn_relu_out        (bn_relu_out), 
    .save_valid         (save_valid),
    .pw_1_bn_relu_done  (pw_1_bn_relu_done)
);

    
//////////////////////////////////////////////////
// BRAM (INPUT, OUTPUT, WEIGHT)

// BRAM_A
blk_mem_gen_0 bram_A (
  .clka                 (clk),          // input wire clka
  .ena                  (ena_0),        // input wire ena
  .wea                  (wea_0),        // input wire [0 : 0] wea
  .addra                (addra_0),      // input wire [8 : 0] addra
  .dina                 (dina_0),       // input wire [3527 : 0] dina
  .clkb                 (clk),          // input wire clkb
  .enb                  (enb_0),        // input wire enb
  .addrb                (addrb_0),      // input wire [8 : 0] addrb
  .doutb                (doutb_0)       // output wire [3527 : 0] doutb
);

// BRAM_B
blk_mem_gen_1 bram_B (
  .clka                 (clk),          // input wire clka
  .ena                  (ena_1),        // input wire ena
  .wea                  (wea_1),        // input wire [0 : 0] wea
  .addra                (addra_1),      // input wire [8 : 0] addra
  .dina                 (dina_1),       // input wire [3527 : 0] dina
  .clkb                 (clk),          // input wire clkb
  .enb                  (enb_1),        // input wire enb
  .addrb                (addrb_1),      // input wire [8 : 0] addrb
  .doutb                (doutb_1)       // output wire [3527 : 0] doutb
);

//////////////////////////////////////////////////////////
// BRAM_W0 (Pointwise_1 Weight)
blk_mem_gen_2 bram_w0 (
  .clka                 (clk),          // input wire clka
  .ena                  (ena_w0),       // input wire ena
  .addra                (addra_w0),     // input wire [14 : 0] addra
  .douta                (douta_w0)      // output wire [16 : 0] douta
);

// BRAM_W1 (Depthwise Weight)
blk_mem_gen_3 bram_w1 (
  .clka                 (clk),          // input wire clka
  .ena                  (ena_w1),       // input wire ena
  .addra                (addra_w1),     // input wire [11 : 0] addra
  .douta                (douta_w1)      // output wire [16 : 0] douta
);


// BRAM_W2 (Pointwise_2 Weight)
blk_mem_gen_4 bram_w2 (
  .clka                 (clk),          // input wire clka
  .ena                  (ena_w2),       // input wire ena
  .addra                (addra_w2),     // input wire [14 : 0] addra
  .douta                (douta_w2)      // output wire [16 : 0] douta
);
///////////////////////////////////////////////////////////////
// bram_PARAMETER 
bram_bias bram_bias_0 (
  .clka                 (clk),              // input wire clka
  .ena                  (ena_bias_0),       // input wire ena
  .addra                (addra_bias_0),     // input wire [8 : 0] addra
  .douta                (douta_bias_0)      // output wire [31 : 0] douta
);

bram_mean bram_mean_0 (
  .clka                 (clk),              // input wire clka
  .ena                  (ena_mean_0),       // input wire ena
  .addra                (addra_mean_0),     // input wire [8 : 0] addra
  .douta                (douta_mean_0)      // output wire [31 : 0] douta
);

bram_std bram_std_0 (
  .clka                 (clk),              // input wire clka
  .ena                  (ena_std_0),        // input wire ena
  .addra                (addra_std_0),      // input wire [8 : 0] addra
  .douta                (douta_std_0)       // output wire [31 : 0] douta
);

bram_weight bram_weight_0 (
  .clka                 (clk),              // input wire clka
  .ena                  (ena_weight_0),     // input wire ena
  .addra                (addra_weight_0),   // input wire [8 : 0] addra
  .douta                (douta_weight_0)    // output wire [31 : 0] douta
);

/////////////////////////////////////////////////////////////

endmodule