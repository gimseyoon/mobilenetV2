`timescale 1ns / 1ps

module tb_mobilenetV2 #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),
    parameter ADDR_WMEM = $clog2(384 * 64),
    parameter ADDR_W1_MEM = $clog2(384 * 9)
);

    reg clk;
    reg rst_n;
    reg start;
    wire result_save_valid;
    wire [13:0] layer_8_result;

    // DUT
    mobilenetV2 #(
        .IO_WIDTH(IO_WIDTH),
        .ROW(ROW),
        .COLUMN(COLUMN),
        .PIXEL(PIXEL),
        .W_WIDTH(W_WIDTH),
        .ADDR_CHANNEL(ADDR_CHANNEL),
        .ADDR_WMEM(ADDR_WMEM)
    ) dut (
        .clk(clk),
        .rst_n(rst_n),
        .start(start),
        .result_save_valid_o(result_save_valid),
        .layer_8_result(layer_8_result)
    );

    // Clock: 300 MHz 예) 주기 3.333...ns -> #1.6666667 토글
    initial clk = 0;
    always #10 clk = ~clk;

    // Reset
    initial begin
        rst_n = 0;
        #25 rst_n = 1;
    end

    // Start
    initial begin
        start = 0;
        #400 start = 1;
        #25  start = 0;
    end

    // ----------------------------
    // 파일로 출력 (skip_valid=1일 때만)
    // ----------------------------
    integer fd;
    integer ch;           // 0..63 (총 64채널)
    integer y, x, idx;
    integer val;          // 32-bit signed for print
    reg result_save_valid_d;     // 에지 검출용

    initial begin
        fd = $fopen("C:/seyoon/mobilenetV2/output.txt", "w");
        if (fd == 0) begin
            $display("ERROR: cannot open output.txt");
            $finish;
        end
        ch = 0;
        result_save_valid_d = 0;
    end

    // skip_valid의 rising-edge에서 한 채널(14x14) 기록
    always @(posedge clk) begin
        result_save_valid_d <= result_save_valid;

        if (rst_n && (result_save_valid && !result_save_valid_d)) begin
            // 채널 헤더(원하면 주석 해제)
            // $fwrite(fd, "CHANNEL %0d\n", ch);

            for (y = 0; y < ROW; y = y + 1) begin
                for (x = 0; x < COLUMN; x = x + 1) begin
                    idx = y*COLUMN + x;
                    // 18-bit slice를 부호확장하여 정수로
                    val = $signed(layer_8_result[IO_WIDTH*(idx+1)-1 -: IO_WIDTH]);
                    // 같은 줄에 공백으로 구분, 줄 끝에 개행
                    if (x == COLUMN-1)
                        $fwrite(fd, "%0d\n", val);
                    else
                        $fwrite(fd, "%0d ",  val);
                end
            end
            // 채널 사이에 빈 줄 하나 (원하지 않으면 주석)
            $fwrite(fd, "\n");

            ch = ch + 1;
            if (ch == 64) begin
                $display("All 64 channels written to output.txt");
                $fclose(fd);
                $finish;
            end
        end
    end

endmodule
