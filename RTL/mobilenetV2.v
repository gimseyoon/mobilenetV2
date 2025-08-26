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
    input clk_in,
    //input clk,
    input rst,
    input start,
    //output result_save_valid_o,
    //output signed [3527:0] result_o,
    output reg [13: 0] layer_8_result,
    output reg all_done
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


    wire clk;
    wire locked;
    wire rst_n_sync;
    

    
/////////////////////////////////////////////////////////////
// FSM
    wire [2:0] state;

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
    
    assign result_save_valid_o = result_save_valid;
    assign result_o = result;
//////////////////////////////////////////////////
// bram_A
    wire                                        ena_0;
    wire                                        wea_0;
    wire  [ADDR_IN-1:0]                         addra_0;
    wire  signed [IO_WIDTH * PIXEL - 1 : 0]     dina_0;
    wire                                        enb_0;
    wire  [ADDR_IN-1:0]                         addrb_0;
    wire  signed [IO_WIDTH * PIXEL - 1 : 0]     doutb_0;
// bram_B
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
    wire                                        ena_biassubam_0;
    wire    [ADDR_PARAM-1 : 0]                  addra_biassubam_0;
    wire    signed      [31:0]                  douta_biassubam_0;
// bram_mean
    wire                                        ena_wdivstd_0;
    wire    [ADDR_PARAM-1 : 0]                  addra_wdivstd_0;
    wire    signed      [31:0]                  douta_wdivstd_0;
    
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
    reg                                         ena_biassubam_0_q;
    reg  [ADDR_PARAM-1 : 0]                     addra_biassubam_0_q;
    reg                                         ena_wdivstd_0_q;
    reg  [ADDR_PARAM-1 : 0]                     addra_wdivstd_0_q;
    
    
    
    
    reg signed [IO_WIDTH -1 : 0] result_1; always@(posedge clk or negedge rst_n_sync) if(!rst_n_sync) result_1 <= 0; else result_1 <= result[IO_WIDTH -1 : 0];
    reg signed [IO_WIDTH -1 : 0] result_2; always@(posedge clk or negedge rst_n_sync) if(!rst_n_sync) result_2 <= 0; else result_2 <= result[2*IO_WIDTH-1 : IO_WIDTH];
    reg signed [IO_WIDTH -1 : 0] result_3; always@(posedge clk or negedge rst_n_sync) if(!rst_n_sync) result_3 <= 0; else result_3 <= result[IO_WIDTH*PIXEL-1-IO_WIDTH -: IO_WIDTH];
    reg signed [IO_WIDTH -1 : 0] result_4; always@(posedge clk or negedge rst_n_sync) if(!rst_n_sync) result_4 <= 0; else result_4 <= result[IO_WIDTH*PIXEL-1 -: IO_WIDTH];
    /*reg all_done;*/                          always@(posedge clk or negedge rst_n_sync) if(!rst_n_sync) all_done <= 0; else all_done <= skip_done;
    
    
////////////////////////////////////////////////////////////
// 1-clk pipeline registers for BRAM/param signals
always @(posedge clk or negedge rst_n_sync) begin
  if (!rst_n_sync) begin
    {ena_0_q,wea_0_q,addra_0_q,dina_0_q,enb_0_q,addrb_0_q,
     ena_1_q,wea_1_q,addra_1_q,dina_1_q,enb_1_q,addrb_1_q,
     ena_w0_q,addra_w0_q, ena_w1_q,addra_w1_q,ena_w2_q,addra_w2_q,
     ena_biassubam_0_q,addra_biassubam_0_q,
     ena_wdivstd_0_q,addra_wdivstd_0_q} <= 0;
  end else begin
    {ena_0_q,wea_0_q,addra_0_q,dina_0_q,enb_0_q,addrb_0_q,
     ena_1_q,wea_1_q,addra_1_q,dina_1_q,enb_1_q,addrb_1_q,
     ena_w0_q,addra_w0_q, ena_w1_q,addra_w1_q,ena_w2_q,addra_w2_q,
     ena_biassubam_0_q,addra_biassubam_0_q,
     ena_wdivstd_0_q,addra_wdivstd_0_q} <=
    {ena_0,wea_0,addra_0,dina_0,enb_0,addrb_0,
     ena_1,wea_1,addra_1,dina_1,enb_1,addrb_1,
     ena_w0,addra_w0, ena_w1,addra_w1,ena_w2,addra_w2,
     ena_biassubam_0,addra_biassubam_0,
     ena_wdivstd_0,addra_wdivstd_0};
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


//////////////////////////////////////////////////
// sk_in_1, sk_in_2    
    assign sk_in_1 = (skip_valid) ? doutb_0 : 0;
    assign sk_in_2 = (skip_valid) ? bn_relu_out : 0;
    

/////////////////////////////////////////////////////// 
// Data muxing for BRAM write data
///////////////////////////////////////////////////////

    assign dina_0 = result_q;
    assign dina_1 = bn_relu_out_q; 
    
//////////////////////////////////////////////////


// ============================================================
// FANOUT SPLIT for state  (no latency)
//    - state_buf : keep
//    - state_row[r] : row ?? ???? (Vivado?? ??? ????????? max_fanout ???)
// ============================================================
(* keep = "true" *) wire [2:0] state_buf = state;
(* max_fanout = 64 *) wire [2:0] state_row [0:ROW-1];

genvar rr;
generate
  for (rr=0; rr<ROW; rr=rr+1) begin : STATE_FANOUT_SPLIT
    assign state_row[rr] = state_buf; // ???? ????? ????(???? ????)
  end
endgenerate

// ============================================================
// decide multiplier input  (row???? select ?? ????/?????? ???)
// ============================================================
integer k;
always @(posedge clk or negedge rst_n_sync) begin
  if(!rst_n_sync) begin
    mul_in     <= 0;
    mul_weight <= 0;
  end else begin
    // ?????? ???(row?? case)
    for (k = 0; k < PIXEL; k = k + 1) begin : MULIN_UPDATE
      case (state_row[k/14])
        PW_1: mul_in[IO_WIDTH*(k+1)-1 -: IO_WIDTH] <= doutb_0[IO_WIDTH*(k+1)-1 -: IO_WIDTH];
        DW  : mul_in[IO_WIDTH*(k+1)-1 -: IO_WIDTH] <= doutb_1[IO_WIDTH*(k+1)-1 -: IO_WIDTH];
        PW_2: mul_in[IO_WIDTH*(k+1)-1 -: IO_WIDTH] <= doutb_1[IO_WIDTH*(k+1)-1 -: IO_WIDTH];
        default: mul_in[IO_WIDTH*(k+1)-1 -: IO_WIDTH] <= 0;
      endcase
    end

    // weight ????(???? ???? ???? select ????)
    case (state_buf)
      PW_1: mul_weight <= douta_w0;
      DW  : mul_weight <= douta_w1;
      PW_2: mul_weight <= douta_w2;
      default: mul_weight <= 0;
    endcase
  end
end


    // start rising-edge ????
    reg start_d;
    always @(posedge clk or negedge rst_n_sync) begin
        if (!rst_n_sync) start_d <= 1'b0;
        else     start_d <= start;
    end
    wire start_rise = start & ~start_d;   // 1clk ???




// Instantiate FSM
reset_sync reset_sync_0(
    .clk(clk),
    .rst(rst),
    .rst_sync(rst_sync),    
    .rst_n_sync(rst_n_sync)
);

//////////////////////////////////////////////////
// Instantiate FSM
FSM FSM_0 (
    .clk                (clk),
    .rst_n              (rst_n_sync),
    .start              (start_rise),
    .pw_1_bn_relu_done  (pw_1_bn_relu_done),
    .dw_bn_relu_done    (dw_bn_relu_done),
    .skip_done          (skip_done),  
    .state              (state)
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
    .ena_biassubam_0    (ena_biassubam_0),
    .addra_biassubam_0  (addra_biassubam_0),
// BRAM mean
    .ena_wdivstd_0      (ena_wdivstd_0),
    .addra_wdivstd_0    (addra_wdivstd_0)
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
    .wdivstd            (douta_wdivstd_0),
    .biassubam          (douta_biassubam_0),
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
bram_B bram_B (
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
bram_W0 bram_w0 (
  .clka                 (clk),          // input wire clka
  .ena                  (ena_w0_q),       // input wire ena
  .addra                (addra_w0_q),     // input wire [14 : 0] addra
  .douta                (douta_w0)      // output wire [16 : 0] douta
);

// BRAM_W1 (Depthwise Weight)
bram_W1 bram_w1 (
  .clka                 (clk),          // input wire clka
  .ena                  (ena_w1_q),       // input wire ena
  .addra                (addra_w1_q),     // input wire [11 : 0] addra
  .douta                (douta_w1)      // output wire [16 : 0] douta
);


// BRAM_W2 (Pointwise_2 Weight)
bram_W2 bram_w2 (
  .clka                 (clk),          // input wire clka
  .ena                  (ena_w2_q),       // input wire ena
  .addra                (addra_w2_q),     // input wire [14 : 0] addra
  .douta                (douta_w2)      // output wire [16 : 0] douta
);
///////////////////////////////////////////////////////////////
// bram_PARAMETER 

bram_bias_sub_am bram_bias_sub_am_0 (
  .clka(clk),    // input wire clka
  .ena(ena_biassubam_0_q),      // input wire ena
  .addra(addra_biassubam_0_q),  // input wire [9 : 0] addra
  .douta(douta_biassubam_0)  // output wire [31 : 0] douta
);

bram_w_div_std bram_w_div_std_0 (
  .clka(clk),    // input wire clka
  .ena(ena_wdivstd_0_q),      // input wire ena
  .addra(addra_wdivstd_0_q),  // input wire [9 : 0] addra
  .douta(douta_wdivstd_0)  // output wire [31 : 0] douta
);
/////////////////////////////////////////////////////////////


clk_wiz_0 clk_100_0 (
    // Clock out ports
    .clk_100(clk),     // output clk_100
    // Status and control signals
    .reset(rst), // input reset
    .locked(locked),       // output locked
   // Clock in ports
    .clk_in1(clk_in)
);      // input clk_in1



/*
ila_1 ila_1_test (
	.clk(clk), // input wire clk


	.probe0(result_save_valid), // input wire [0:0]  probe0  
	.probe1(start), // input wire [0:0]  probe1 
	.probe2(result_1), // input wire [17:0]  probe2 
	.probe3(result_4), // input wire [17:0]  probe3 
	.probe4(state), // input wire [2:0]  probe4 
	.probe5(all_done), // input wire [0:0]  probe5 
	.probe6(start_d), // input wire [0:0]  probe6 
	.probe7(start_rise), // input wire [0:0]  probe7 
	.probe8(pw_1_bn_relu_done), // input wire [0:0]  probe8 
	.probe9(dw_bn_relu_done), // input wire [0:0]  probe9 
	.probe10(skip_done), // input wire [0:0]  probe10 
	.probe11(skip_valid), // input wire [0:0]  probe11 
	.probe12(save_valid), // input wire [0:0]  probe12 
	.probe13(pw_1_valid), // input wire [0:0]  probe13 
	.probe14(pw_1_done), // input wire [0:0]  probe14 
	.probe15(dw_valid), // input wire [0:0]  probe15 
	.probe16(dw_done), // input wire [0:0]  probe16 
	.probe17(pw_2_valid), // input wire [0:0]  probe17 
	.probe18(pw_2_done), // input wire [0:0]  probe18 
	.probe19(clk), // input wire [0:0]  probe19
	.probe20(rst), // input wire [0:0]  probe20 
	.probe21(rst_sync), // input wire [0:0]  probe21 
	.probe22(rst_n_sync) // input wire [0:0]  probe22
);
*/

endmodule