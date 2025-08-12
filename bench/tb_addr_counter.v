`timescale 1ns / 1ps

module tb_addr_counter;

    // Parameters
    parameter IO_WIDTH       = 18;
    parameter ROW            = 14;
    parameter COLUMN         = 14;
    parameter PIXEL          = ROW * COLUMN;
    parameter W_WIDTH        = 17;
    parameter ADDR_CHANNEL   = $clog2(384);        // 9
    parameter ADDR_WMEM      = $clog2(384 * 64);   // 15
    parameter ADDR_W1_MEM    = $clog2(384 * 9);    // 12

    // DUT I/O
    reg                         clk;
    reg                         rst_n;
    reg [2:0]                   state;
    reg [14:0]                  glbl_cnt;
    reg [ADDR_CHANNEL-1:0]      acc_cnt;
    reg                         save_valid;
    reg                         pw_1_valid;
    reg                         enb_0;
    reg                         enb_1;

    wire [ADDR_CHANNEL-1:0]     channel_num;
    wire [ADDR_CHANNEL-1:0]     addra_0;
    wire [ADDR_CHANNEL-1:0]     addrb_0;
    wire [ADDR_CHANNEL-1:0]     addra_1;
    wire [ADDR_CHANNEL-1:0]     addrb_1;
    wire [ADDR_WMEM-1:0]        addra_w0;
    wire [ADDR_W1_MEM-1:0]      addra_w1;
    wire [ADDR_WMEM-1:0]        addra_w2;
    wire [ADDR_CHANNEL-1:0]     addra_bias_0;
    wire [ADDR_CHANNEL-1:0]     addra_mean_0;
    wire [ADDR_CHANNEL-1:0]     addra_std_0;
    wire [ADDR_CHANNEL-1:0]     addra_weight_0;

    // Instantiate DUT
    addr_counter #(
        .IO_WIDTH(IO_WIDTH),
        .ROW(ROW),
        .COLUMN(COLUMN),
        .PIXEL(PIXEL),
        .W_WIDTH(W_WIDTH),
        .ADDR_CHANNEL(ADDR_CHANNEL),
        .ADDR_WMEM(ADDR_WMEM),
        .ADDR_W1_MEM(ADDR_W1_MEM)
    ) uut (
        .clk(clk),
        .rst_n(rst_n),
        .state(state),
        .glbl_cnt(glbl_cnt),
        .acc_cnt(acc_cnt),
        .save_valid(save_valid),
        .pw_1_valid(pw_1_valid),
        .enb_0(enb_0),
        .enb_1(enb_1),
        .channel_num(channel_num),
        .addra_0(addra_0),
        .addrb_0(addrb_0),
        .addra_1(addra_1),
        .addrb_1(addrb_1),
        .addra_w0(addra_w0),
        .addra_w1(addra_w1),
        .addra_w2(addra_w2),
        .addra_bias_0(addra_bias_0),
        .addra_mean_0(addra_mean_0),
        .addra_std_0(addra_std_0),
        .addra_weight_0(addra_weight_0)
    );

    // Clock generation
    initial clk = 0;
    always #5 clk = ~clk; // 100MHz

    // Test logic
    initial begin
        $display("Starting TB...");
        rst_n = 0;
        state = 3'b000;
        glbl_cnt = 0;
        acc_cnt = 0;
        save_valid = 0;
        pw_1_valid = 0;
        enb_0 = 0;
        enb_1 = 0;

        #20 rst_n = 1;
        #10;

        // Set to DW mode
        state = 3'b011; // DW
        enb_1 = 1;

        repeat (40) begin
            #10;
            glbl_cnt = glbl_cnt + 1;
            $display("[%0t ns] glbl_cnt=%0d, addrb_1=%0d", $time, glbl_cnt, addrb_1);
        end

        $display("TB finished.");
        $stop;
    end

<<<<<<< HEAD
endmodule
=======
endmodule
>>>>>>> cd5d375ebbe5b2d89c1741d5fc379a7068d5594e
