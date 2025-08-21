`timescale 1ns / 1ps

module tb_mobilenetV2#(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),
    parameter ADDR_WMEM = $clog2(384 * 64),
    parameter ADDR_W1_MEM = $clog2(384 * 9)
);

  reg clk = 1'b0;
  reg rst = 1'b1;
  reg start = 1'b0;

  // 모니터링: ReLU 결과
  wire signed [3527:0] result;
  wire [13:0] layer_8_result;
  wire result_save_valid;
    // Clock generation
    always #5 clk = ~clk; // 100MHz (10ns period)

    // Stimulus
    initial begin
        #200;
        rst = 0;
        #20;
        rst = 1;
        #20;

        start = 1;
        #10;
        start = 0;

        // 충분한 처리 시간 확보
        #900000;

        #200;
        $finish;
    end



  // -------------------------
  // DUT
  // -------------------------
  mobilenetV2 DUT (
    .clk(clk),
    .rst_n(rst),
    .start(start),
    .result_save_valid_o(result_save_valid),
    .result_o(result),
    .layer_8_result(layer_8_result)
  );
  

integer fd;

initial begin
    fd = $fopen("result.txt", "w");
    if (fd == 0) begin
        $display("파일 열기 실패");
        $finish;
    end
end

// 1) 우선 동작 확인용: 매 클럭마다 시간과 sum_out을 찍어보기
//    (잘 나오면 아래의 if(in_valid) 조건 버전으로 바꾸세요)
always @(posedge clk) begin
    if(result_save_valid) begin
        $fdisplay(fd, "%0t %04h", $time, result);
    end
end

// 2) 시뮬 끝나기 직전에 닫기 (같은 초기 블록에서 $finish보다 먼저!)
initial begin
    #500000000;             // 사용자가 말한 시뮬 종료 시점
    $fclose(fd);       // 반드시 먼저 닫고
    $finish;           // 그 다음 종료
end
endmodule

    

