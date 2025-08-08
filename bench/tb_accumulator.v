`timescale 1ns / 1ps

module tb_accumulator;

    localparam IO_WIDTH = 18;
    localparam ROW = 14;
    localparam COLUMN = 14;
    localparam PIXEL = ROW * COLUMN;
    localparam W_WIDTH = 17;
    localparam ADDR_CHANNEL  = $clog2(384);
    localparam ADDR_WMEM = $clog2(384 * 64);
    localparam ADDR_W1_MEM = $clog2(384 * 9);

    // DUT signals
    reg clk;
    reg rst_n;
    reg [2:0] state;
    reg signed [IO_WIDTH * PIXEL - 1:0] mul_out;

    wire [ADDR_CHANNEL -1 :0] acc_cnt;
    wire pw_1_valid;
    wire pw_1_done;
    wire dw_valid;
    wire dw_done;
    wire pw_2_valid;
    wire pw_2_done;
    wire [ADDR_CHANNEL-1 :0] channel_num;
    wire signed [IO_WIDTH * PIXEL - 1 : 0] acc_out;

    // DUT instantiation
    accumulator #(
        .IO_WIDTH(IO_WIDTH),
        .ROW(ROW),
        .COLUMN(COLUMN),
        .PIXEL(PIXEL),
        .W_WIDTH(W_WIDTH),
        .ADDR_CHANNEL(ADDR_CHANNEL),
        .ADDR_WMEM(ADDR_WMEM),
        .ADDR_W1_MEM(ADDR_W1_MEM)
    ) dut (
        .clk(clk),
        .rst_n(rst_n),
        .state(state),
        .mul_out(mul_out),
        .acc_cnt(acc_cnt),
        .pw_1_valid(pw_1_valid),
        .pw_1_done(pw_1_done),
        .dw_valid(dw_valid),
        .dw_done(dw_done),
        .pw_2_valid(pw_2_valid),
        .pw_2_done(pw_2_done),
        .channel_num(channel_num),
        .acc_out(acc_out)
    );

    // Clock generation
    initial clk = 0;
    always #5 clk = ~clk;

    // File logging
    integer outfile;
    integer i, r, c;
    reg signed [IO_WIDTH-1:0] acc_matrix [0:PIXEL-1];

    // Simulation control
initial begin
    // 1. 초기화
    rst_n = 0;
    state = 3'b000;
    mul_out = 0;
    outfile = $fopen("C:/seyoon/acc_out_result.txt", "w");
    if (!outfile) begin
        $display("ERROR: Cannot open output file.");
        $finish;
    end

    // 2. Reset & 상태 설정
    #20;
    rst_n = 1;
    #15;

    state = 3'b011; // DW 상태 진입
    #50;

    // 3. mul_out 모든 값 1로 설정
    for (i = 0; i < PIXEL; i = i + 1) begin
        mul_out[i*IO_WIDTH +: IO_WIDTH] = 1;
    end
end

// acc_out 값을 95ns부터 1클럭마다 출력
initial begin
    #95;  // 정확히 95ns 후에 시작

    for (i = 0; i < 9; i = i + 1) begin
        @(posedge clk);
        
        for (r = 0; r < PIXEL; r = r + 1)
            acc_matrix[r] = acc_out[r*IO_WIDTH +: IO_WIDTH];

        $fdisplay(outfile, "=== acc_out at clk %0d (time: %0t ns) ===", i, $time);
        for (r = 0; r < ROW; r = r + 1) begin
            for (c = 0; c < COLUMN; c = c + 1) begin
                $fwrite(outfile, "%0d ", acc_matrix[r * COLUMN + c]);
            end
            $fwrite(outfile, "\n");
        end
        $fwrite(outfile, "\n");
    end

    $fclose(outfile);
    $finish;
end

endmodule