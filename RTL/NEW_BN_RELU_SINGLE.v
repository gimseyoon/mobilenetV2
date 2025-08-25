`timescale 1ns / 1ps

module NEW_BN_RELU_SINGLE #(
    parameter IO_WIDTH = 18
)(
    input                       clk,
    input                       rst_n,
    input                       bn_en,

    input       [31:0]          biassubam,
    input       [31:0]          wdivstd,
    
    input  signed [IO_WIDTH-1:0] acc_in,
    output signed [IO_WIDTH-1:0] bn_out,
    output                       valid_out
);

/////////////////////////////////////////////////////// 
// Wires between FP IP blocks
///////////////////////////////////////////////////////
wire signed [31:0] data_f2f_out;
wire               valid_f2f_out;

wire signed [31:0] data_sub_out;
wire               valid_sub_out;

wire signed [31:0] data_div_out;
wire               valid_div_out;

wire signed [31:0] data_mul_out;
wire               valid_mul_out;

wire signed [31:0] data_add_out;
wire               valid_add_out;

wire signed [IO_WIDTH-1:0] data_f2fx_out;
wire                       valid_f2fx_out;

/////////////////////////////////////////////////////// 
// Outputs
///////////////////////////////////////////////////////
assign bn_out    = valid_f2fx_out ? data_f2fx_out : 0;
assign valid_out = valid_f2fx_out;

/////////////////////////////////////////////////////// 
// Fixed -> Float (acc_in)
///////////////////////////////////////////////////////
fixed_to_float fixed_to_float_0 (
    .aclk                 (clk),
    .aresetn              (rst_n),
    .s_axis_a_tvalid      (bn_en),
    .s_axis_a_tdata       (acc_in),
    .m_axis_result_tvalid (valid_f2f_out),
    .m_axis_result_tdata  (data_f2f_out)
);


/////////////////////////////////////////////////////// 
// Multiply: * wdivstd
///////////////////////////////////////////////////////
Multiply Multiply_0 (
    .aclk                 (clk),
    .aresetn              (rst_n),
    .s_axis_a_tvalid      (valid_f2f_out),
    .s_axis_a_tdata       (data_f2f_out),
    .s_axis_b_tvalid      (valid_f2f_out),
    .s_axis_b_tdata       (wdivstd),
    .m_axis_result_tvalid (valid_mul_out),
    .m_axis_result_tdata  (data_mul_out)
);

/////////////////////////////////////////////////////// 
// Add: + biassubam
///////////////////////////////////////////////////////
add add_0 (
    .aclk                 (clk),
    .aresetn              (rst_n),
    .s_axis_a_tvalid      (valid_mul_out),
    .s_axis_a_tdata       (data_mul_out),
    .s_axis_b_tvalid      (valid_mul_out),
    .s_axis_b_tdata       (biassubam),
    .m_axis_result_tvalid (valid_add_out),
    .m_axis_result_tdata  (data_add_out)
);

/////////////////////////////////////////////////////// 
// Float -> Fixed (1.2.15 into IO_WIDTH)
///////////////////////////////////////////////////////
float_to_fixed float_to_fixed_0 (
    .aclk                 (clk),
    .aresetn              (rst_n),
    .s_axis_a_tvalid      (valid_add_out),
    .s_axis_a_tdata       (data_add_out),
    .m_axis_result_tvalid (valid_f2fx_out),
    .m_axis_result_tdata  (data_f2fx_out)
);

endmodule