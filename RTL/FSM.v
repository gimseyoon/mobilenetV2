`timescale 1ns / 1ps

module FSM(
    input               clk,
    input               rst_n,
    input               start,
    input               new_start,
    input               pw_1_bn_relu_done,
    input               dw_bn_relu_done,
    input               skip_done,     // layer_8_done
    output reg  [2:0]  state,
    output reg  [1:0]  layer_state
);

/////////////////////////////////////////////////////// 
// State encoding
///////////////////////////////////////////////////////
localparam IDLE     = 3'b000,
           PW_1     = 3'b001,
           PW_1_RST = 3'b010,
           DW       = 3'b011,
           DW_RST   = 3'b100,
           PW_2     = 3'b101,
           PW_2_RST = 3'b110,
           EXPORT   = 3'b111;

localparam READY    = 2'b00,
           LAYER_8  = 2'b01,
           LAYER_9  = 2'b10,
           LAYER_10 = 2'b11;
           
/////////////////////////////////////////////////////// 
// Next-state register
///////////////////////////////////////////////////////
reg [2:0] next_state;
reg [1:0] next_layer_state;
/////////////////////////////////////////////////////// 
// State register
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state <= IDLE;
        layer_state <= READY;
    end
    else begin
        state <= next_state;
        layer_state <= next_layer_state;
    end
end

/////////////////////////////////////////////////////// 
// Next-state combinational logic
///////////////////////////////////////////////////////
always @(*) begin
    case (state)
        IDLE:     if (start || new_start)   next_state = PW_1; else next_state = next_state;
        PW_1:     if (pw_1_bn_relu_done)    next_state = PW_1_RST; else next_state = next_state;
        PW_1_RST:                           next_state = DW;
        DW:       if (dw_bn_relu_done)      next_state = DW_RST; else next_state = next_state;
        DW_RST:                             next_state = PW_2;
        PW_2:     if (skip_done)            next_state = IDLE; else next_state = next_state;
        /*
        PW_2_RST:                           next_state = EXPORT;
        EXPORT:   if (export_done)          next_state = IDLE;
        */
        default: next_state = IDLE;
    endcase
end

always @(*) begin
    case (layer_state)
        READY:       if (start || new_start)    next_layer_state = LAYER_8; else next_layer_state = next_layer_state; 
        LAYER_8:     if (new_start)             next_layer_state = LAYER_9; else next_layer_state = next_layer_state; 
        LAYER_9:     if (new_start)             next_layer_state = LAYER_10; else next_layer_state = next_layer_state; 
        LAYER_10:    if (new_start)             next_layer_state = READY; else next_layer_state = next_layer_state; 

        default: next_layer_state = READY;
    endcase
end



endmodule