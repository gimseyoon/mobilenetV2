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
