module BN_RELU #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,                 // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),          // 8 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64)          // 15 (for 64*384 = 24576)
)(
    input                                       clk,
    input                                       rst_n,
    input                                       bn_relu_valid,
    input                                       mean,
    input                                       weight,
    input                                       bias,
    input                                       std,
    input signed [IO_WIDTH * PIXEL - 1 : 0]     acc_out,        // [3920-1 : 0], 3920 bit
    output signed [IO_WIDTH * PIXEL - 1 : 0]    bn_relu_out,    // [3920-1 : 0], 3920 bit
    output                                      save_valid
    );
    
//////////////////////////////////////////////////

    // Fixed to Float
    wire signed [IO_WIDTH-1 : 0] data_f2f_in  = acc_out;
    wire signed [31:0] data_f2f_out;
    wire               valid_f2f_out;
    // Subtract
    wire signed [31:0] data_sub_out;
    wire               valid_sub_out;
    // Divide
    wire signed [31:0] data_div_out;
    wire               valid_div_out;
    // Multiply
    wire signed [31:0] data_mul_out;
    wire               valid_mul_out;
    // Add
    wire signed [31:0] data_add_out;
    wire               valid_add_out;
    // Float to Fixed
    wire signed [IO_WIDTH-1 : 0] data_f2fx_out;
    wire               valid_f2fx_out;
    
    
//////////////////////////////////////////////////

    assign bn_relu_out = (valid_f2fx_out) ? data_f2fx_out[15:0] : 0;
    assign save_valid = valid_f2fx_out;
    
    
    
//////////////////////////////////////////////////
// IP Connections


// Fixed to Float
fixed_to_float fixed_to_float_0 (
  .aclk(clk),
  .aresetn(rst_n),
  .s_axis_a_tvalid(valid_f2f_in),
  .s_axis_a_tdata(data_f2f_in),
  .m_axis_result_tvalid(valid_f2f_out),
  .m_axis_result_tdata(data_f2f_out)
);

// Subtract: x - mean
subtract subtract_0 (
  .aclk(clk),
  .aresetn(rst_n),
  .s_axis_a_tvalid(valid_f2f_out),
  .s_axis_a_tdata(data_f2f_out),
  .s_axis_b_tvalid(valid_f2f_out),
  .s_axis_b_tdata(mean),
  .m_axis_result_tvalid(valid_sub_out),
  .m_axis_result_tdata(data_sub_out)
);

// Divide: (x - mean) / std
divide divide_0 (
  .aclk(clk),
  .aresetn(rst_n),
  .s_axis_a_tvalid(valid_sub_out),
  .s_axis_a_tdata(data_sub_out),
  .s_axis_b_tvalid(valid_sub_out),
  .s_axis_b_tdata(std),
  .m_axis_result_tvalid(valid_div_out),
  .m_axis_result_tdata(data_div_out)
);

// Multiply: * weight
Multiply Multiply_0 (
  .aclk(clk),
  .aresetn(rst_n),
  .s_axis_a_tvalid(valid_div_out),
  .s_axis_a_tdata(data_div_out),
  .s_axis_b_tvalid(valid_div_out),
  .s_axis_b_tdata(weight),
  .m_axis_result_tvalid(valid_mul_out),
  .m_axis_result_tdata(data_mul_out)
);

// Add: + bias
add add_0 (
  .aclk(clk),
  .aresetn(rst_n),
  .s_axis_a_tvalid(valid_mul_out),
  .s_axis_a_tdata(data_mul_out),
  .s_axis_b_tvalid(valid_mul_out),
  .s_axis_b_tdata(bias),
  .m_axis_result_tvalid(valid_add_out),
  .m_axis_result_tdata(data_add_out)
);

// Float to Fixed
float_to_fixed float_to_fixed_0 (
  .aclk(clk),
  .aresetn(rst_n),
  .s_axis_a_tvalid(valid_add_out),
  .s_axis_a_tdata(data_add_out),
  .m_axis_result_tvalid(valid_f2fx_out),
  .m_axis_result_tdata(data_f2fx_out)
);


endmodule
