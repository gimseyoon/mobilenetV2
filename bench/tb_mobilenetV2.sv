`timescale 1ns / 1ps

module tb_mobilenetV2#(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),
    parameter ADDR_WMEM = $clog2(384 * 64),
    parameter ADDR_W1_MEM = $clog2(384 * 9)
);

  reg clk = 1'b0;
  reg rst_n = 1'b1;
  reg start = 1'b0;

  wire signed [3527:0] result;
  wire result_save_valid;
  
  reg [17:0] result_d [0:195];
  reg result_save_valid_d;
  
  
  
  //wire all_done;
  //wire [1:0] layer_8_result;
  
    // Clock generation
    always #5 clk = ~clk; // 100MHz (10ns period)

    // Stimulus
    initial begin
        #200;
        rst_n = 0;
        #20;
        rst_n = 1;
        #20;

        #1000;
        start = 1;
        #10;
        start = 0;


    end




  // -------------------------
  // DUT
  // -------------------------
  mobilenetV2 DUT (
    .clk(clk),
    .rst_n(rst_n),
    .start(start),
    .result_save_valid_o(result_save_valid),
    .result_o(result),
    //.all_done(all_done),
    .layer_8_result(layer_8_result)
  );

integer p;
always@(posedge clk or negedge rst_n) begin
    if(!rst_n) begin
        for (p = 0; p < PIXEL; p = p + 1) begin
            result_d[p] <= 0;
        end
        result_save_valid_d <= 0;
    end
    else begin
    
    for (p = 0; p < PIXEL; p = p + 1) begin
        result_d[p] <= result[p*IO_WIDTH +: IO_WIDTH];
    end  
        result_save_valid_d <= result_save_valid;
        
    end
end




integer fd;

initial begin
    fd = $fopen("result.txt", "w");
    if (fd == 0) begin
        $display("파일 열기 실패");
        $finish;
    end
end

// ---- 블록 카운트 ----
integer block_cnt = 0;

// result_d 샘플 & valid 딜레이는 기존 그대로 사용한다고 가정
// result_save_valid_d 가 1일 때마다 196개 덤프
always @(posedge clk) begin
    if (result_save_valid_d) begin
        integer i;
        reg [19:0] hex20;
        for (i = 0; i < PIXEL; i = i + 1) begin
            hex20 = { {2{result_d[i][17]}}, result_d[i] }; // 18b -> 20b 부호확장
            $fwrite(fd, "%05h", hex20);
            if (i != PIXEL-1) $fwrite(fd, " ");
        end
        $fwrite(fd, "\n\n");

        // ---- 블록 수 증가 ----
        block_cnt <= block_cnt + 1;

        // ---- 64번째 쓴 "바로 그 사이클"에 종료 ----
        if (block_cnt == 63) begin
            // block_cnt가 이번에 64가 되므로(0부터 시작) 바로 종료
            $fflush(fd);   // 버퍼 강제 플러시
            $fclose(fd);
            $display("[%0t] 모든 64 블록 기록 완료", $time);
            $finish;
        end
    end
end



endmodule

    

