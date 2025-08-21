// multiplier
// - 18b × [196] mul_in 을 17b weight와 곱해 18b로 round+shift
// - weight는 행/소그룹 단위로 레지스터 복제해 fanout과 net delay 완화
`timescale 1ns / 1ps

module multiplier #(
    parameter IO_WIDTH        = 18,
    parameter ROW             = 14,
    parameter COLUMN          = 14,
    parameter PIXEL           = ROW * COLUMN,   // 196
    parameter W_WIDTH         = 17,
    parameter integer R_SHIFT = 16,
    // 한 행을 SUB개 소그룹으로 나눠 weight 레지스터 복제 (1,2,4…)
    parameter integer SUB     = 2
)(
    input                               clk,
    input                               rst_n,
    input  signed [IO_WIDTH*PIXEL-1:0]  mul_in,
    input  signed [W_WIDTH-1:0]         mul_weight,
    output signed [IO_WIDTH*PIXEL-1:0]  mul_out
);

///////////////////////////////////////////////////////
// 입력 1clk 파이프라인
///////////////////////////////////////////////////////
reg signed [IO_WIDTH*PIXEL-1:0] mul_in_q;
always @(posedge clk or negedge rst_n)
    if (!rst_n) mul_in_q <= {IO_WIDTH*PIXEL{1'b0}};
    else        mul_in_q <= mul_in;

///////////////////////////////////////////////////////
// weight 복제 레지스터 (평면 1D 배열: 인덱스 = 행*SUB + 그룹)
///////////////////////////////////////////////////////
(* DONT_TOUCH = "true" *)
reg signed [W_WIDTH-1:0] w_dup [0:ROW*SUB-1];

integer rr, gg;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        for (rr=0; rr<ROW*SUB; rr=rr+1) w_dup[rr] <= {W_WIDTH{1'b0}};
    end else begin
        // 동일 weight를 행×소그룹 수만큼 복제(레지스터화)
        for (rr=0; rr<ROW; rr=rr+1)
            for (gg=0; gg<SUB; gg=gg+1)
                w_dup[rr*SUB + gg] <= mul_weight;
    end
end

///////////////////////////////////////////////////////
// 곱 결과 라운드/시프트 후 레지스터
///////////////////////////////////////////////////////
function signed [17:0] round_shift_signed18;
    input signed [34:0] x;
    reg sign; reg [34:0] mag; reg [35:0] mag_rnd; reg [17:0] q_u;
begin
    sign    = x[34];
    mag     = sign ? (~x + 35'd1) : x;
    mag_rnd = {1'b0, mag} + (36'd1 << (R_SHIFT-1));
    q_u     = mag_rnd[35:R_SHIFT];
    round_shift_signed18 = sign ? -$signed(q_u) : $signed(q_u);
end
endfunction

reg  signed [IO_WIDTH-1:0] mul_out_reg [0:PIXEL-1];
wire signed [34:0]         mul_out_w   [0:PIXEL-1];

integer k;
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        for (k=0; k<PIXEL; k=k+1) mul_out_reg[k] <= {IO_WIDTH{1'b0}};
    end else begin
        for (k=0; k<PIXEL; k=k+1) mul_out_reg[k] <= round_shift_signed18(mul_out_w[k]);
    end
end

genvar m;
generate
    for (m=0; m<PIXEL; m=m+1) begin : PACK
        assign mul_out[IO_WIDTH*(m+1)-1 : IO_WIDTH*m] = mul_out_reg[m];
    end
endgenerate

///////////////////////////////////////////////////////
// DSP 인스턴스: B는 복제 레지스터 w_dup 사용(교차-참조 없음)
///////////////////////////////////////////////////////
genvar ry, cx;
generate
    for (ry=0; ry<ROW; ry=ry+1) begin : ROWS
        for (cx=0; cx<COLUMN; cx=cx+1) begin : COLS
            localparam integer IDX  = ry*COLUMN + cx;
            // 열을 SUB개 그룹으로 균등 분할해 선택(정수 상수 → 정적 인덱스)
            localparam integer GIDX = (cx*SUB)/COLUMN; // 0..SUB-1

            (* use_dsp = "yes" *)
            mult_gen_0 u_mul (
                .CLK (clk),
                .A   ($signed(mul_in_q[IO_WIDTH*(IDX+1)-1 : IO_WIDTH*IDX])),
                .B   (w_dup[ry*SUB + GIDX]),
                .P   (mul_out_w[IDX])
            );
        end
    end
endgenerate

endmodule
