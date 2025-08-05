`timescale 1ns / 1ps

module mobilenetV2#(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),         // CHANNEL = 384
    parameter ADDR_WMEM = $clog2(384 * 64)       // 15 (for 24576)
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
    wire signed [IO_WIDTH * PIXEL - 1 : 0]   acc_out;
    wire save_valid;
    wire [8:0] channel_num;
    wire pw_1_done;
    
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

//////////////////////////////////////////////////

    assign mul_weight = douta_w0;
    assign mul_in = doutb_0;
    


//////////////////////////////////////////////////
// Instantiate FSM
FSM FSM_0(
    .clk(clk),
    .rst_n(rst_n),
    .start(start),
    .pw_1_done(pw_1_done),
    .pw_1_bn_relu_done(pw_1_bn_relu_done),
    .dw_done(dw_done),
    .dw_bn_relu_done(dw_bn_relu_done),
    .pw_2_done(pw_2_done),
    .pw_2_bn_done(pw_2_bn_done),
    .layer_done(layer_done),  
    .bram_select(bram_select),
    .state(state)
);
    

//////////////////////////////////////////////////
// Instantiate addr_counter

glbl_ctrl glbl_ctrl_0 (
    .clk          (clk),
    .rst_n        (rst_n),
    .state        (state),
    .save_valid   (save_valid),
    .channel_num  (channel_num),
    .bram_select  (bram_select),
    .acc_out      (acc_out),

    // BRAM 0
    .ena_0        (ena_0),
    .wea_0        (wea_0),
    .addra_0      (addra_0),
    .dina_0       (dina_0),
    .enb_0        (enb_0),
    .addrb_0      (addrb_0),

    // BRAM 1
    .ena_1        (ena_1),
    .wea_1        (wea_1),
    .addra_1      (addra_1),
    .dina_1       (dina_1),    
    .enb_1        (enb_1),
    .addrb_1      (addrb_1),

    // BRAM W
    .ena_w0       (ena_w0),
    .addra_w0     (addra_w0)
);



    
//////////////////////////////////////////////////
// Instantiate accumulator
accumulator #(
    .IO_WIDTH(IO_WIDTH),
    .ROW(ROW),
    .COLUMN(COLUMN),
    .PIXEL(PIXEL),
    .W_WIDTH(W_WIDTH)
) accumulator_0 (
    .clk(clk),
    .rst_n(rst_n),
    .state(state),
    .mul_out(mul_out),
    .save_valid(save_valid),
    .channel_num(channel_num),
    .pw_1_done(pw_1_done),
    .acc_out(acc_out)
);
    
//////////////////////////////////////////////////
// Instantiate multiplier
multiplier #(
    .IO_WIDTH(IO_WIDTH),
    .ROW(ROW),
    .COLUMN(COLUMN),
    .PIXEL(PIXEL),
    .W_WIDTH(W_WIDTH)
) multiplier_0 (
    .clk(clk),
    .rst_n(rst_n),
    .mul_in(mul_in),
    .mul_weight(mul_weight),
    .mul_out(mul_out)
);
    
//////////////////////////////////////////////////
// Instantiate BRAM

blk_mem_gen_0 IO_bram_A (
  .clka(clk),    // input wire clka
  .ena(ena_0),      // input wire ena
  .wea(wea_0),      // input wire [0 : 0] wea
  .addra(addra_0),  // input wire [8 : 0] addra
  .dina(dina_0),    // input wire [3527 : 0] dina
  .clkb(clk),    // input wire clkb
  .enb(enb_0),      // input wire enb
  .addrb(addrb_0),  // input wire [8 : 0] addrb
  .doutb(doutb_0)  // output wire [3527 : 0] doutb
);


blk_mem_gen_1 IO_bram_B (
  .clka(clk),    // input wire clka
  .ena(ena_1),      // input wire ena
  .wea(wea_1),      // input wire [0 : 0] wea
  .addra(addra_1),  // input wire [8 : 0] addra
  .dina(dina_1),    // input wire [3527 : 0] dina
  .clkb(clk),    // input wire clkb
  .enb(enb_1),      // input wire enb
  .addrb(addrb_1),  // input wire [8 : 0] addrb
  .doutb(doutb_1)  // output wire [3527 : 0] doutb
);


blk_mem_gen_2 W_bram_w0 (
  .clka(clk),    // input wire clka
  .ena(ena_w0),      // input wire ena
  .addra(addra_w0),  // input wire [14 : 0] addra
  .douta(douta_w0)  // output wire [16 : 0] douta
);

endmodule