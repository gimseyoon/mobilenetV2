`timescale 1ns / 1ps
module reset_sync #(
    parameter integer STAGES = 2  // 권장: 2 이상
)(
    input  wire clk,
    input  wire rst,          // 보드에서 들어오는 비동기 리셋 (Active-HIGH, posedge)
    output wire rst_n_sync,   // 동기 해제된 리셋 (Active-LOW)  -> 모듈의 negedge 포트에 연결
    output wire rst_sync      // 동기 해제된 리셋 (Active-HIGH) -> 필요시 사용
);
    // STAGES 최소값 체크 (시뮬/합성 경고 방지)
    initial begin
        if (STAGES < 2) $error("reset_sync: STAGES must be >= 2");
    end

    // 비동기 리셋용 메타/동기화 레지스터
    // Xilinx 권장 속성: ASYNC_REG (배치 시 메타 FF로 인식)
    (* ASYNC_REG = "TRUE" *) reg [STAGES-1:0] shreg;

    // 비동기 assert (posedge rst), 동기 de-assert (clk)
    always @(posedge clk or posedge rst) begin
        if (rst) begin
            // 리셋이 들어오면 즉시 LOW로 끌어내려서 (Active-LOW 출력 기준)
            // 모든 하위 모듈이 비동기적으로 리셋되도록 함
            shreg <= {STAGES{1'b0}};
        end else begin
            // 클럭 도메인 안에서 한 스테이지씩 1로 채워 넣으며
            // de-assert를 동기화
            shreg <= {shreg[STAGES-2:0], 1'b1};
        end
    end

    // 최종 출력: Active-LOW / Active-HIGH 모두 제공
    assign rst_n_sync = shreg[STAGES-1];      // 모듈의 negedge 리셋 포트에 바로 연결
    assign rst_sync   = ~rst_n_sync;          // 필요시 Active-HIGH 형태

endmodule
