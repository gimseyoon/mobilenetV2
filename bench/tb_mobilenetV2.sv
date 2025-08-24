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
  reg rst = 1'b0;
  reg start = 1'b0;

  wire signed [3527:0] result;
  wire result_save_valid;
  
  reg [17:0] result_d [0:195];
  reg result_save_valid_d;
  
  
  
  //wire all_done;
  //wire [1:0] layer_8_result;
  
    // Clock generation
    always #2.5 clk = ~clk; // 200MHz (5ns period)

    // Stimulus
    initial begin
        #200;
        rst = 1;
        #10;
        rst = 0;
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
    .rst(rst),
    .start(start),
    .result_save_valid_o(result_save_valid),
    .result_o(result),
    .all_done(all_done),
    .layer_8_result(layer_8_result)
  );

integer p;
always@(posedge clk or negedge rst) begin
    if(rst) begin
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
    fd = $fopen("result_all_layer_onlydata.txt", "w");
    if (fd == 0) begin
        $display("파일 열기 실패");
        $finish;
    end
end

// ---- 블록 카운트 ----
integer block_cnt = 0;


always @(posedge clk) begin
    if (result_save_valid_d) begin
        integer i;
        integer this_blk;        // 이번 블록 인덱스(0~191)
        integer layer_idx;       // 8/9/10 (참고용)
        integer chan_idx;        // 0~63   (참고용)
        reg [19:0] hex20;

        this_blk  = block_cnt;
        layer_idx = 8 + (this_blk / 64);   // 8,9,10
        chan_idx  = this_blk % 64;         // 0..63

        // ===== 헤더 =====
        // 요구사항: [CHANNEL N]
        // (원하면 아래처럼 레이어/채널로도 출력 가능)
        // $fwrite(fd, "[LAYER %0d] [CHANNEL %0d]\n", layer_idx, chan_idx);

        // ===== 값 196개(18b → 부호확장 20b, 공백 구분) =====
        for (i = 0; i < PIXEL; i = i + 1) begin
            hex20 = { {2{result_d[i][17]}}, result_d[i] };
            $fwrite(fd, "%05h", hex20);
            if (i != PIXEL-1) $fwrite(fd, " ");
        end
        $fwrite(fd, "\n");   // 줄바꿈 (헤더 뒤 1줄 + 값 뒤 1줄 = 총 2줄)

        // ===== 카운터 갱신 및 종료 처리 =====
        block_cnt = block_cnt + 1;
        if (this_blk == 191) begin           // 0~191 = 192번째를 막 쓴 사이클
            $fflush(fd);
            $fclose(fd);
            $display("[%0t] 192 blocks (LAYER8~10) saved.", $time);
            $finish;
        end
    end
end



endmodule

    

