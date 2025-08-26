`timescale 1ns / 1ps
module SK #(
    parameter IO_WIDTH = 18,
    parameter ROW      = 14,
    parameter COLUMN   = 14,
    parameter PIXEL    = ROW * COLUMN            // 196
)(
    input                                   clk,
    input                                   rst_n,
    input                                   skip_valid,        
    input       signed [IO_WIDTH*PIXEL-1:0] in_1,
    input       signed [IO_WIDTH*PIXEL-1:0] in_2,
    output      signed [IO_WIDTH*PIXEL-1:0] result,                   
    output reg                              result_save_valid,         
    output reg                              skip_done               
);

    // -------------------------
    reg        busy;
    reg  [7:0] idx;                        // 0..195
    reg        skip_valid_d;
    wire       start = skip_valid & ~skip_valid_d;

    reg  [PIXEL-1:0] wr_sel;             
    wire             last_write = wr_sel[PIXEL-1];

    reg  [6:0] sk_cnt;                 
    reg        skip_done_pending;        


    // -------------------------
    reg  signed [IO_WIDTH*PIXEL-1:0] in_1_q, in_2_q;
    wire signed [IO_WIDTH-1:0] in1_sel = in_1_q[IO_WIDTH*idx +: IO_WIDTH];
    wire signed [IO_WIDTH-1:0] in2_sel = in_2_q[IO_WIDTH*idx +: IO_WIDTH];
    wire signed [IO_WIDTH-1:0] sum_w   = in1_sel + in2_sel;  

 
 
    // -------------------------
    reg  signed [IO_WIDTH-1:0] result_reg [0:PIXEL-1];


    // -------------------------
    wire [IO_WIDTH*PIXEL-1:0] packed_result;
    reg  [IO_WIDTH*PIXEL-1:0] result_q;
    assign result = result_q;

    genvar g;
    generate
        for (g = 0; g < PIXEL; g = g + 1) begin : PACK_RESULT
            assign packed_result[IO_WIDTH*(g+1)-1 : IO_WIDTH*g] = result_reg[g];
        end
    endgenerate


    reg save_pending;

    integer k;

    // -------------------------

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            busy               <= 1'b0;
            idx                <= 8'd0;
            wr_sel             <= {PIXEL{1'b0}};
            skip_valid_d       <= 1'b0;
            result_save_valid  <= 1'b0;
            skip_done          <= 1'b0;
            sk_cnt             <= 7'd0;
            skip_done_pending  <= 1'b0;
            in_1_q             <= {IO_WIDTH*PIXEL{1'b0}};
            in_2_q             <= {IO_WIDTH*PIXEL{1'b0}};
            save_pending       <= 1'b0;
            result_q           <= {IO_WIDTH*PIXEL{1'b0}};
            for (k = 0; k < PIXEL; k = k + 1) result_reg[k] <= {IO_WIDTH{1'b0}};
        end else begin
            skip_valid_d      <= skip_valid;
            result_save_valid <= 1'b0;
            skip_done         <= 1'b0;


            if (start && !busy) begin
                in_1_q <= in_1;
                in_2_q <= in_2;
                busy   <= 1'b1;
                idx    <= 8'd0;
                wr_sel <= {{(PIXEL-1){1'b0}}, 1'b1};  // LSB=1
            end

            else if (busy) begin

                for (k = 0; k < PIXEL; k = k + 1)
                    if (wr_sel[k]) result_reg[k] <= sum_w;

   
                if (last_write) begin
                    busy              <= 1'b0;
                    wr_sel            <= {PIXEL{1'b0}};
                    save_pending      <= 1'b1;         
                    skip_done_pending <= (sk_cnt == 7'd63); 


                    if (sk_cnt == 7'd63) sk_cnt <= 7'd0;
                    else                  sk_cnt <= sk_cnt + 1'b1;
                end else begin

                    wr_sel <= {wr_sel[PIXEL-2:0], 1'b0};
                    idx    <= idx + 1'b1;
                end
            end


            if (save_pending) begin
                result_q          <= packed_result;   
                result_save_valid <= 1'b1;            
                skip_done         <= skip_done_pending;
                skip_done_pending <= 1'b0;
                save_pending      <= 1'b0;
            end
        end
    end
endmodule