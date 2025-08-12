<<<<<<< HEAD
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
=======
`timescale 1ns/1ps

module tb_accumulator;

  localparam IO_WIDTH = 18;
  localparam ROW = 14;
  localparam COLUMN = 14;
  localparam PIXEL = ROW*COLUMN;

  reg  clk;
  reg  rst_n;
  reg  [2:0] state;
  reg  signed [IO_WIDTH*PIXEL-1:0] mul_out;

  wire bn_en;
  wire [8:0] acc_cnt;
  wire pw_1_valid, pw_1_done;
  wire dw_valid,  dw_done;
  wire pw_2_valid, pw_2_done;
  wire [8:0] channel_num;
  wire signed [IO_WIDTH*PIXEL-1:0] acc_out;

  // DUT
  accumulator #(
    .IO_WIDTH(IO_WIDTH), .ROW(ROW), .COLUMN(COLUMN)
  ) dut (
    .clk(clk), .rst_n(rst_n),
    .state(state),
    .mul_out(mul_out),
    .bn_en(bn_en),
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

  // clock
  initial clk = 0;
  always #5 clk = ~clk;

  // util: set mul_out pixels
  task set_all_pixels(input signed [IO_WIDTH-1:0] v);
    integer i;
    begin
      for (i=0; i<PIXEL; i=i+1)
        mul_out[i*IO_WIDTH +: IO_WIDTH] = v;
    end
  endtask

  integer fd, y, x;


  task dump_frame_acc_out_reg(input integer step_id);
    integer idx;
    begin
      $fwrite(fd, "===== STEP %0d (acc_cnt=%0d) =====\n", step_id, acc_cnt);
      for (y=0; y<ROW; y=y+1) begin
        for (x=0; x<COLUMN; x=x+1) begin
          idx = y*COLUMN + x;

          $fwrite(fd, "%0d", $signed(dut.acc_out_reg[idx]));
          if (x != COLUMN-1) $fwrite(fd, " ");
        end
        $fwrite(fd, "\n");
      end
      $fwrite(fd, "\n");
    end
  endtask


  task dump_frame_acc_out_bus;
    integer idx;
    begin
      $fwrite(fd, "===== FINAL (dw_valid=1) acc_out bus =====\n");
      for (y=0; y<ROW; y=y+1) begin
        for (x=0; x<COLUMN; x=x+1) begin
          idx = y*COLUMN + x;
          $fwrite(fd, "%0d", $signed(acc_out[(idx+1)*IO_WIDTH-1 -: IO_WIDTH]));
          if (x != COLUMN-1) $fwrite(fd, " ");
        end
        $fwrite(fd, "\n");
      end
      $fwrite(fd, "\n");
    end
  endtask

  integer step;

    initial begin
    rst_n = 0;
    state = 3'b000; // IDLE
    mul_out = '0;

    fd = $fopen("C:/seyoon/mobilenetV2/value/result/reference_result/layer_8/test.txt", "w");
    if (fd == 0) begin
      $display("ERROR: cannot open output file.");
      $finish;
    end

    // reset
    repeat (10) @(posedge clk);
    rst_n = 1;


    set_all_pixels(18'sd1);


    @(posedge clk);
    state = 3'b011; // DW

    wait (dut.dw_en == 1'b1);


    wait (acc_cnt == 9'd0);
    @(posedge clk);
    wait (acc_cnt == 9'd1);


    for (step = 1; step <= 9; step = step + 1) begin
      wait (acc_cnt == step[8:0]);
      @(posedge clk);              
      dump_frame_acc_out_reg(step);
    end


    $fclose(fd);
    $display("[TB] first set (acc_cnt 1..9) saved.");
    $finish;
  end

endmodule
>>>>>>> cd5d375ebbe5b2d89c1741d5fc379a7068d5594e
