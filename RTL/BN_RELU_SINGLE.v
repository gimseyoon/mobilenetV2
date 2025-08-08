module BN_RELU_SINGLE #(
    parameter IO_WIDTH = 18
)(
    input clk,
    input rst_n,
    input pw_1_valid,
    input bn_en,

    input [31:0] mean,
    input [31:0] weight,
    input [31:0] bias,
    input [31:0] std,

    input  signed [IO_WIDTH-1:0] acc_in,
    output signed [IO_WIDTH-1:0] bn_out,
    output valid_out
);

/////////////////////////////////////////////////////////////////////////////

    wire signed [31:0] data_f2f_out;
    wire valid_f2f_out;

    wire signed [31:0] data_sub_out;
    wire valid_sub_out;

    wire signed [31:0] data_div_out;
    wire valid_div_out;

    wire signed [31:0] data_mul_out;
    wire valid_mul_out;

    wire signed [31:0] data_add_out;
    wire valid_add_out;

    wire signed [IO_WIDTH-1:0] data_f2fx_out;
    wire valid_f2fx_out;

/////////////////////////////////////////////////////////////////////////////

    assign bn_out     = (valid_f2fx_out) ? data_f2fx_out[15:0] : 0;
    assign valid_out  = valid_f2fx_out;
    
/////////////////////////////////////////////////////////////////////////////
    
    // Fixed to Float
    fixed_to_float fixed_to_float_0 (
        .aclk(clk),
        .aresetn(rst_n),
        .s_axis_a_tvalid(bn_en),
        .s_axis_a_tdata(acc_in),
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

    // Divide
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

    // Multiply
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

    // Add
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

/////////////////////////////////////////////////////////////////////////////
    
    
    
endmodule
