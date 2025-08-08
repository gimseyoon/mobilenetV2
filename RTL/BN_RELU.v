module BN_RELU #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,                 // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),          // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64)          // 15 (for 64*384 = 24576)
)(
    input                                       clk,
    input                                       rst_n,
    input                                       pw_1_valid,
    input                                       bn_en,
    input                             [31:0]    mean,
    input                             [31:0]    weight,
    input                             [31:0]    bias,
    input                             [31:0]    std,
    input  signed [IO_WIDTH * PIXEL - 1 : 0]    acc_out,        // [3920-1 : 0], 3920 bit
    output signed [IO_WIDTH * PIXEL - 1 : 0]    bn_relu_out,    // [3920-1 : 0], 3920 bit
    output                                      save_valid,
    output reg                                 pw_1_bn_relu_done
    );
    
//////////////////////////////////////////////////
    reg                    [5:0] bn_cnt;            // 0~48
    reg      [ADDR_CHANNEL-1 :0] bn_channel_num;
    reg  signed [IO_WIDTH-1 : 0] acc_out_reg       [0 : PIXEL-1];
    wire        [IO_WIDTH-1 : 0] acc_selected      [3:0];
    wire         [IO_WIDTH-1 :0] bn_single_out     [3:0];
    reg          [IO_WIDTH-1 :0] bn_relu_out_array [PIXEL-1 :0];




    
///////////////////////////////////////////////////////////////////////
// assign acc_out_array[i] = acc_out[20(i+1) -1 :20*i]

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        bn_cnt <= 0;
    end
    else begin
        if (bn_cnt == 48) begin 
            bn_cnt <= 0;
        end
        else if(bn_en) begin
            bn_cnt <= bn_cnt + 1;
        end
    end
end


///////////////////////////////////////////////////////////////////////
// assign acc_out_array[i] = acc_out[20(i+1) -1 :20*i]

genvar i;
generate
    for (i = 0; i < PIXEL; i = i + 1) begin : REG_ARRAY_ASSIGN
        always @(posedge clk or negedge rst_n) begin
            if (!rst_n)
                acc_out_reg[i] <= 0;
            else
                if(pw_1_valid) begin
                    acc_out_reg[i] <= acc_out[IO_WIDTH*(i+1)-1 : IO_WIDTH*i];
                end
        end
    end
endgenerate
    
    
    
///////////////////////////////////////////////////////////////////////
// aacc_selected <= acc_out_reg
genvar s;
generate
    for (s = 0; s < 4; s = s + 1) begin : SELECT_ACC
        assign acc_selected[s] = acc_out_reg[bn_cnt*4 + s];
    end
endgenerate



///////////////////////////////////////////////////////////////////////
// assign bn_relu_out[3920-1 :0] = bn_relu_out_array [20-1 :0][196-1 :0]

    genvar k;
    generate
        for (k = 0; k < PIXEL; k = k + 1) begin : PACK_OUTPUT
            assign bn_relu_out [(IO_WIDTH*k+1)-1 : IO_WIDTH*k] = (save_valid) ? bn_relu_out_array[k] : 0;
        end
    endgenerate
    
    
    
///////////////////////////////////////////////////////////////////////
// channel_num

    always@(posedge clk or negedge rst_n) begin
        if(!rst_n) begin
            bn_channel_num <= 0;
        end
        else begin
            if (save_valid) begin bn_channel_num <= (bn_channel_num == 383) ? 0 : bn_channel_num + 1; end
        end
    end


///////////////////////////////////////////////////////////////////////
// bn_relu_out_array

    integer q;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (q = 0; q < PIXEL; q = q + 1)
                bn_relu_out_array[q] <= 0;
        end 
        else begin
            if (save_valid) begin
                bn_relu_out_array[bn_cnt * 4 + 0] <= bn_single_out[0];
                bn_relu_out_array[bn_cnt * 4 + 1] <= bn_single_out[1];
                bn_relu_out_array[bn_cnt * 4 + 2] <= bn_single_out[2];
                bn_relu_out_array[bn_cnt * 4 + 3] <= bn_single_out[3];
            end 
        end 
    end

///////////////////////////////////////////////////////////////////////
// pw_1_bn_relu_done

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (q = 0; q < PIXEL; q = q + 1)
                pw_1_bn_relu_done <= 0;
        end 
        else begin
            if ( (bn_cnt == 6'd48) && (bn_channel_num == 9'd383) ) begin   
                pw_1_bn_relu_done <= 1;
            end
            else begin
                pw_1_bn_relu_done <= 0;
            end
        end //else
    end //always
    
    
///////////////////////////////////////////////////////////////////////
    genvar p;
    generate
        for (p = 0; p < 4; p = p + 1) begin : BN_UNIT
            BN_RELU_SINGLE #(
                .IO_WIDTH(IO_WIDTH)
            ) BN_RELU_SINGLE_0 (
                .clk(clk),
                .rst_n(rst_n),
                .pw_1_valid(pw_1_valid),
                .bn_en(bn_en),
                .mean(mean),
                .weight(weight),
                .bias(bias),
                .std(std),
                .acc_in(acc_selected[i]),
                .bn_out(bn_single_out[i]),
                .valid_out(save_valid)
            );
        end
    endgenerate
    
    
///////////////////////////////////////////////////////////////////////
endmodule
