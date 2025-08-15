`timescale 1ns / 1ps

module FSM(
    input               clk,
    input               rst_n,
    input               start,
    input               pw_1_done,
    input               pw_1_bn_relu_done,
    input               dw_done,
    input               dw_bn_relu_done,
    input               pw_2_done,
    input               pw_2_bn_done,
    input               layer_done,     // PW_2 done
    output reg  [2:0]   state
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

/////////////////////////////////////////////////////// 
// Next-state register
///////////////////////////////////////////////////////
reg [2:0] ns;

/////////////////////////////////////////////////////// 
// State register
///////////////////////////////////////////////////////
always @(posedge clk or negedge rst_n) begin
    if (!rst_n)
        state <= IDLE;
    else
        state <= ns;
end

/////////////////////////////////////////////////////// 
// Next-state combinational logic
///////////////////////////////////////////////////////
always @(*) begin
    ns = state;
    case (state)
        IDLE:     if (start)            ns = PW_1;
        PW_1:     if (pw_1_bn_relu_done)ns = PW_1_RST;
        PW_1_RST:                       ns = DW;
        DW:       if (dw_bn_relu_done)  ns = DW_RST;
        DW_RST:                         ns = PW_2;
        PW_2:     if (pw_2_bn_done)     ns = PW_2_RST;
        PW_2_RST:                       ns = EXPORT;
        EXPORT:   if (layer_done)       ns = IDLE;
        default:                        ns = IDLE;
    endcase
end

endmodule
