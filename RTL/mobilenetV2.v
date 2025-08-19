`timescale 1ns / 1ps

module mobilenetV2 #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter INPUT_CHANNEL = 64,
    parameter ADDR_PARAM = 10,
    parameter ADDR_IN = $clog2(64),
    parameter ADDR_CHANNEL  = $clog2(384),        // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64),       // 15 (for 64*384 = 24576)
    parameter ADDR_W1_MEM = $clog2(384 * 9)       // 12 (for 9*384 = 3456)
)(
    input clk,
    input rst_n,
    input start,
    //output reg result_save_valid_o,
    output reg [13: 0] layer_8_result

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
              
    localparam READY    = 2'b00,
               LAYER_8  = 2'b01,
               LAYER_9  = 2'b10,
               LAYER_10 = 2'b11;
               
/////////////////////////////////////////////////////////////

    reg new_start;
    wire rst_n_sync;
    
/////////////////////////////////////////////////////////////
// FSM
    wire [2:0] state;
    wire[1:0] layer_state;

//////////////////////////////////////////////////
// multiplier
    reg signed [IO_WIDTH * PIXEL - 1 : 0]   mul_in;   // [3528-1 : 0], 3528 bit
    reg signed [W_WIDTH - 1 : 0]            mul_weight;        // [17-1:0], 17 bit
    wire signed [IO_WIDTH * PIXEL - 1 : 0]  mul_out; // [3528-1 : 0], 3528 bit

//////////////////////////////////////////////////
//accumultaor
    wire signed [IO_WIDTH * PIXEL - 1 : 0] acc_out;
    wire pw_1_valid;   
    wire pw_1_done;
    wire dw_valid;
    wire dw_done;
    wire pw_2_valid;
    wire pw_2_done;

//////////////////////////////////////////////////
// BN_RELU
    wire signed [IO_WIDTH * PIXEL - 1 : 0] bn_relu_out;
    reg  signed [IO_WIDTH * PIXEL - 1 : 0] bn_relu_out_q; 
    wire pw_1_bn_relu_done;
    wire dw_bn_relu_done;
    wire pw_2_bn_done;
    wire save_valid;
    wire skip_valid;
    
    always @(posedge clk) begin
        bn_relu_out_q <= bn_relu_out;
    end
//////////////////////////////////////////////////
// SK
    wire signed [IO_WIDTH * PIXEL - 1 : 0] sk_in_1;
    wire signed [IO_WIDTH * PIXEL - 1 : 0] sk_in_2;
    wire signed [IO_WIDTH * PIXEL - 1 : 0] result;
    reg  signed [IO_WIDTH * PIXEL - 1 : 0] result_q;
    wire result_save_valid;
    wire skip_done;
    
    always @(posedge clk) begin
        result_q <= result;
    end
//////////////////////////////////////////////////
// bram_0
    wire                                        ena_0;
    wire                                        wea_0;
    wire  [ADDR_IN-1:0]                         addra_0;
    wire  signed [IO_WIDTH * PIXEL - 1 : 0]     dina_0;
    wire                                        enb_0;
    wire  [ADDR_IN-1:0]                         addrb_0;
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
    wire    [ADDR_PARAM-1 : 0]                  addra_bias_0;
    wire    signed      [31:0]                  douta_bias_0;
// bram_mean
    wire                                        ena_mean_0;
    wire    [ADDR_PARAM-1 : 0]                  addra_mean_0;
    wire    signed      [31:0]                  douta_mean_0;
// bram_std
    wire                                        ena_std_0;
    wire    [ADDR_PARAM-1 : 0]                  addra_std_0;
    wire    signed      [31:0]                  douta_std_0;
// bram_weight
    wire                                        ena_weight_0;
    wire    [ADDR_PARAM-1 : 0]                  addra_weight_0;
    wire    signed      [31:0]                  douta_weight_0;
    
    //////////////////////////////////////////////////
    // bram_q
    reg                                         ena_0_q;
    reg                                         wea_0_q;
    reg  [ADDR_IN-1:0]                          addra_0_q;
    reg  signed [IO_WIDTH * PIXEL - 1 : 0]      dina_0_q;
    reg                                         enb_0_q;
    reg  [ADDR_IN-1:0]                          addrb_0_q;
    reg                                         ena_1_q;
    reg                                         wea_1_q;
    reg  [ADDR_CHANNEL-1:0]                     addra_1_q;
    reg  signed [IO_WIDTH * PIXEL - 1 : 0]      dina_1_q;
    reg                                         enb_1_q;
    reg  [ADDR_CHANNEL-1:0]                     addrb_1_q;
    reg                                         ena_w0_q;
    reg  [ADDR_WMEM-1:0]                        addra_w0_q; 
    reg                                         ena_w1_q;
    reg  [ADDR_W1_MEM-1:0]                      addra_w1_q; 
    reg                                         ena_w2_q;
    reg  [ADDR_WMEM-1:0]                        addra_w2_q; 
    reg                                         ena_bias_0_q;
    reg  [ADDR_PARAM-1 : 0]                     addra_bias_0_q;
    reg                                         ena_mean_0_q;
    reg  [ADDR_PARAM-1 : 0]                     addra_mean_0_q;
    reg                                         ena_std_0_q;
    reg  [ADDR_PARAM-1 : 0]                     addra_std_0_q;
    reg                                         ena_weight_0_q;
    reg  [ADDR_PARAM-1 : 0]                     addra_weight_0_q;
    
////////////////////////////////////////////////////////////
// 1-clk pipeline registers for BRAM/param signals
always @(posedge clk or negedge rst_n) begin
  if (!rst_n) begin
    {ena_0_q,wea_0_q,addra_0_q,dina_0_q,enb_0_q,addrb_0_q,
     ena_1_q,wea_1_q,addra_1_q,dina_1_q,enb_1_q,addrb_1_q,
     ena_w0_q,addra_w0_q, ena_w1_q,addra_w1_q,ena_w2_q,addra_w2_q,
     ena_bias_0_q,addra_bias_0_q,
     ena_mean_0_q,addra_mean_0_q,
     ena_std_0_q,addra_std_0_q,
     ena_weight_0_q,addra_weight_0_q} <= 0;
  end else begin
    {ena_0_q,wea_0_q,addra_0_q,dina_0_q,enb_0_q,addrb_0_q,
     ena_1_q,wea_1_q,addra_1_q,dina_1_q,enb_1_q,addrb_1_q,
     ena_w0_q,addra_w0_q, ena_w1_q,addra_w1_q,ena_w2_q,addra_w2_q,
     ena_bias_0_q,addra_bias_0_q,
     ena_mean_0_q,addra_mean_0_q,
     ena_std_0_q,addra_std_0_q,
     ena_weight_0_q,addra_weight_0_q} <=
    {ena_0,wea_0,addra_0,dina_0,enb_0,addrb_0,
     ena_1,wea_1,addra_1,dina_1,enb_1,addrb_1,
     ena_w0,addra_w0, ena_w1,addra_w1,ena_w2,addra_w2,
     ena_bias_0,addra_bias_0,
     ena_mean_0,addra_mean_0,
     ena_std_0,addra_std_0,
     ena_weight_0,addra_weight_0};
  end
end



/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

    integer j;
    always@(posedge clk) begin
        for (j = 0; j < 14; j = j + 1) begin
            layer_8_result[j] <= result[(j * 18 * 14) + 17];
        end  
    end
    /*
    always@(posedge clk) begin
            result_save_valid_o <= result_save_valid;
    end
    */
//////////////////////////////////////////////////
// new_start ( layer (9, 10) )
/*
    always@(posedge clk or negedge rst_n) begin
        if(!rst_n) begin
            new_start <= 0;
        end
        else begin
            if(skip_done) begin
                new_start <= 1;
            end
            else begin
                new_start <= 0;
            end
        end // rst_n
    end //always
 */ 
//////////////////////////////////////////////////
// sk_in_1, sk_in_2    
    assign sk_in_1 = (skip_valid) ? doutb_0 : 0;
    assign sk_in_2 = (skip_valid) ? bn_relu_out : 0;
    
//////////////////////////////////////////////////
// decide multiplier input
always @(posedge clk or negedge rst_n) begin
    if(!rst_n) begin
        mul_in <= 0;
    end
    else begin
        case (state)
          PW_1: begin mul_in <= doutb_0; mul_weight <= douta_w0; end
          DW  : begin mul_in <= doutb_1; mul_weight <= douta_w1; end
          PW_2: begin mul_in <= doutb_1; mul_weight <= douta_w2; end
          default: begin mul_in <= 0; mul_weight <= 0; end
        endcase 
    end

end 

/////////////////////////////////////////////////////// 
// Data muxing for BRAM write data
///////////////////////////////////////////////////////

    assign dina_0 = result_q;
    assign dina_1 = bn_relu_out_q; 
    
//////////////////////////////////////////////////


// Instantiate FSM
reset_sync reset_sync_0(
    .clk(clk),
    .rst_n_async(rst_n),     // 보드에서 들어온 비동기 리셋(Active-Low)
    .rst_n_sync(rst_n_sync)
);

//////////////////////////////////////////////////
// Instantiate FSM
FSM FSM_0 (
    .clk                (clk),
    .rst_n              (rst_n_sync),
    .start              (start),
    .new_start          (new_start),
    .pw_1_bn_relu_done  (pw_1_bn_relu_done),
    .dw_bn_relu_done    (dw_bn_relu_done),
    .skip_done          (skip_done),  
    .state              (state),
    .layer_state        (layer_state)
);
    

//////////////////////////////////////////////////
// Instantiate addr_counter

glbl_ctrl glbl_ctrl_0 (
    .clk                (clk),
    .rst_n              (rst_n_sync),
    .state              (state),
    .save_valid         (save_valid),
    .skip_valid         (skip_valid),
    .acc_out            (acc_out),
    .pw_1_done          (pw_1_done),
    .result_save_valid  (result_save_valid),
    
// BRAM A
    .ena_0              (ena_0),
    .wea_0              (wea_0),
    .addra_0            (addra_0),
    .enb_0              (enb_0),
    .addrb_0            (addrb_0),
// BRAM B
    .ena_1              (ena_1),
    .wea_1              (wea_1),
    .addra_1            (addra_1),  
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
    .rst_n              (rst_n_sync),
    .state              (state),
    .mul_out            (mul_out),
    .bn_en              (bn_en),
    .pw_1_valid         (pw_1_valid),
    .pw_1_done          (pw_1_done),
    .dw_valid           (dw_valid),
    .dw_done            (dw_done),
    .pw_2_valid         (pw_2_valid),
    .pw_2_done          (pw_2_done),
    .acc_out            (acc_out)
);
    
//////////////////////////////////////////////////
// Instantiate multiplier
multiplier multiplier_0 (
    .clk                (clk),
    .rst_n              (rst_n_sync),
    .pw_1_done          (pw_1_done),
    .dw_done            (dw_done),
    .pw_2_done          (pw_2_done),
    .mul_in             (mul_in),     
    .mul_weight         (mul_weight),
    .mul_out            (mul_out)
);

//////////////////////////////////////////////////
// Instantiate BN_RELU
BN_RELU BN_RELU_0 (
    .clk                (clk),
    .rst_n              (rst_n_sync),
    .state              (state),
    .pw_1_valid         (pw_1_valid),
    .dw_valid           (dw_valid),
    .pw_2_valid         (pw_2_valid),
    .bn_en              (bn_en),
    .mean               (douta_mean_0),
    .weight             (douta_weight_0),
    .bias               (douta_bias_0),
    .std                (douta_std_0),
    .acc_out            (acc_out), 
    .bn_relu_out        (bn_relu_out), 
    .save_valid         (save_valid),
    .skip_valid         (skip_valid),
    .pw_1_bn_relu_done  (pw_1_bn_relu_done),
    .dw_bn_relu_done    (dw_bn_relu_done),
    .pw_2_bn_done        (pw_2_bn_done)
);

//////////////////////////////////////////////////
// Instantiate SKIP_CONNECTION
SK SK_0(
    .clk(clk),
    .rst_n(rst_n_sync),
    .skip_valid(skip_valid),
    .in_1(sk_in_1),
    .in_2(sk_in_2),
    .result(result),
    .result_save_valid(result_save_valid),
    .skip_done(skip_done)
    );
    
//////////////////////////////////////////////////
// BRAM (INPUT, OUTPUT, WEIGHT)

// BRAM_A

bram_A bram_A (
  .clka                 (clk),    // input wire clka
  .ena                  (ena_0_q),      // input wire ena
  .wea                  (wea_0_q),      // input wire [0 : 0] wea
  .addra                (addra_0_q),  // input wire [5 : 0] addra
  .dina                 (dina_0),    // input wire [3527 : 0] dina
  .clkb                 (clk),    // input wire clkb
  .enb                  (enb_0_q),      // input wire enb
  .addrb                (addrb_0_q),  // input wire [5 : 0] addrb
  .doutb                (doutb_0)  // output wire [3527 : 0] doutb
);
// BRAM_B
blk_mem_gen_1 bram_B (
  .clka                 (clk),          // input wire clka
  .ena                  (ena_1_q),        // input wire ena
  .wea                  (wea_1_q),        // input wire [0 : 0] wea
  .addra                (addra_1_q),      // input wire [8 : 0] addra
  .dina                 (dina_1),       // input wire [3527 : 0] dina
  .clkb                 (clk),          // input wire clkb
  .enb                  (enb_1_q),        // input wire enb
  .addrb                (addrb_1_q),      // input wire [8 : 0] addrb
  .doutb                (doutb_1)       // output wire [3527 : 0] doutb
);

//////////////////////////////////////////////////////////
// BRAM_W0 (Pointwise_1 Weight)
blk_mem_gen_2 bram_w0 (
  .clka                 (clk),          // input wire clka
  .ena                  (ena_w0_q),       // input wire ena
  .addra                (addra_w0_q),     // input wire [14 : 0] addra
  .douta                (douta_w0)      // output wire [16 : 0] douta
);

// BRAM_W1 (Depthwise Weight)
blk_mem_gen_3 bram_w1 (
  .clka                 (clk),          // input wire clka
  .ena                  (ena_w1_q),       // input wire ena
  .addra                (addra_w1_q),     // input wire [11 : 0] addra
  .douta                (douta_w1)      // output wire [16 : 0] douta
);


// BRAM_W2 (Pointwise_2 Weight)
blk_mem_gen_4 bram_w2 (
  .clka                 (clk),          // input wire clka
  .ena                  (ena_w2_q),       // input wire ena
  .addra                (addra_w2_q),     // input wire [14 : 0] addra
  .douta                (douta_w2)      // output wire [16 : 0] douta
);
///////////////////////////////////////////////////////////////
// bram_PARAMETER 
bram_bias bram_bias_0 (
  .clka(clk),    // input wire clka
  .ena(ena_bias_0_q),      // input wire ena
  .addra(addra_bias_0_q),  // input wire [9 : 0] addra
  .douta(douta_bias_0)  // output wire [31 : 0] douta
);

bram_mean bram_mean_0 (
  .clka(clk),    // input wire clka
  .ena(ena_mean_0_q),      // input wire ena
  .addra(addra_mean_0_q),  // input wire [9 : 0] addra
  .douta(douta_mean_0)  // output wire [31 : 0] douta
);

bram_std bram_std_0 (
  .clka(clk),    // input wire clka
  .ena(ena_std_0_q),      // input wire ena
  .addra(addra_std_0_q),  // input wire [9 : 0] addra
  .douta(douta_std_0)  // output wire [31 : 0] douta
);

bram_weight bram_weight_0 (
  .clka(clk),    // input wire clka
  .ena(ena_weight_0_q),      // input wire ena
  .addra(addra_weight_0_q),  // input wire [9 : 0] addra
  .douta(douta_weight_0)  // output wire [31 : 0] douta
);
/////////////////////////////////////////////////////////////

endmodule