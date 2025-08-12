`timescale 1ns/1ps

module tb_bn_relu;

  // ------------------------------------------------------------
  // DUT 파라미터
  // ------------------------------------------------------------
  localparam IO_WIDTH    = 18;
  localparam ROW         = 14;
  localparam COLUMN      = 14;
  localparam PIXEL       = ROW * COLUMN;          // 196
  localparam ADDR_CH_AW  = $clog2(384);
  localparam ADDR_WMEM   = $clog2(384*64);

  // ------------------------------------------------------------
  // TB 신호
  // ------------------------------------------------------------
  reg                          clk;
  reg                          rst_n;
  reg                          pw_1_valid;
  reg                          bn_en;

  // BN 파라미터(IEEE754 float)
  // mean=0, std=1.0, weight=1.0, bias=0
  reg  [31:0]                  mean;
  reg  [31:0]                  std;
  reg  [31:0]                  weight;
  reg  [31:0]                  bias;

  reg  signed [IO_WIDTH*PIXEL-1:0] acc_out;       // 3920-bit 입력 버스
  wire signed [IO_WIDTH*PIXEL-1:0] bn_relu_out;   // 3920-bit 출력 버스
  wire                            save_valid;
  wire                            pw_1_bn_relu_done;

  // ------------------------------------------------------------
  // DUT 인스턴스
  // ------------------------------------------------------------
  BN_RELU #(
    .IO_WIDTH(IO_WIDTH),
    .ROW(ROW),
    .COLUMN(COLUMN)
  ) dut (
    .clk(clk),
    .rst_n(rst_n),
    .pw_1_valid(pw_1_valid),
    .bn_en(bn_en),
    .mean(mean),
    .weight(weight),
    .bias(bias),
    .std(std),
    .acc_out(acc_out),
    .bn_relu_out(bn_relu_out),
    .save_valid(save_valid),
    .pw_1_bn_relu_done(pw_1_bn_relu_done)
  );

  // ------------------------------------------------------------
  // 클럭/리셋
  // ------------------------------------------------------------
  initial clk = 1'b0;
  always #5 clk = ~clk;  // 100MHz

  initial begin
    rst_n = 1'b0;
    pw_1_valid = 1'b0;
    bn_en = 1'b0;

    mean   = 32'h00000000; // 0.0
    std    = 32'h3F800000; // 1.0
    weight = 32'h3F800000; // 1.0
    bias   = 32'h00000000; // 0.0

    acc_out = 0;

    repeat (10) @(posedge clk);
    rst_n = 1'b1;

    // ----------------------------------------------------------
    // 입력 프레임 한 채널 채우기
    // acc_out[ IO_WIDTH*(i+1)-1 : IO_WIDTH*i ] = 어떤 패턴
    // 여기선 간단히 i를 18비트로 집어넣음(부호 고려해 양수범위)
    // ----------------------------------------------------------
    fill_acc_with_pattern();

    // acc_out_reg에 적재: pw_1_valid 1클럭
    @(posedge clk);
    pw_1_valid <= 1'b1;
    @(posedge clk);
    pw_1_valid <= 1'b0;

    // bn_cnt 굴리기: 파이프라인이 있으므로 지속적으로 enable
    bn_en <= 1'b1;

    // 충분히 돈 뒤 관찰
    repeat (2000) @(posedge clk);

    $display("[TB] done.");
    $finish;
  end

  // ------------------------------------------------------------
  // 유틸: acc_out 버스 채우기
  // ------------------------------------------------------------
  task fill_acc_with_pattern;
    integer i;
    reg signed [IO_WIDTH-1:0] val;
    begin
      for (i = 0; i < PIXEL; i = i + 1) begin
        // 예시 패턴: 0.0 ~ 0.5 근처 값이 되도록 (Q1.17 가정)
        // IO_WIDTH=18, signed. 0.048 * 2^17 ? 6291 정도
        // 다양한 값:  (i%32)*500 정도로 스윕
        val = (i%32) * 500;   // 부호 양수, 약한 증가
        set_pixel(i, val);
      end
    end
  endtask

  // acc_out의 i번째 픽셀 슬라이스에 값 쓰기
  // i번째 픽셀에 값 쓰기 (변수 인덱스 가능)
task set_pixel(input integer idx, input signed [IO_WIDTH-1:0] v);
  begin
    acc_out[idx*IO_WIDTH +: IO_WIDTH] = v;
  end
endtask

  // ------------------------------------------------------------
  // 간단 모니터
  // ------------------------------------------------------------
  // 첫 8개 픽셀만 샘플 출력(필요 시 늘려서 확인)
  integer j;
  always @(posedge clk) begin
    if (save_valid) begin
      $display("save_valid=1 : bn_relu_out[0..7] = ");
      for (j=0; j<8; j=j+1) begin
        $write("%0d ", $signed(bn_relu_out[IO_WIDTH*(j+1)-1 -: IO_WIDTH]));
      end
      $write("\n");
    end
  end

endmodule
