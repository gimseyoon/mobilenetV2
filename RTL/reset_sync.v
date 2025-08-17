module reset_sync #(
    parameter STAGES = 2
)(
    input  wire clk,
    input  wire rst_n_async,     // 보드에서 들어온 비동기 리셋(Active-Low)
    output wire rst_n_sync
);
    reg [STAGES-1:0] shreg;
    always @(posedge clk or negedge rst_n_async) begin
        if (!rst_n_async)
            shreg <= {STAGES{1'b0}};           // 비동기 assert
        else
            shreg <= {shreg[STAGES-2:0], 1'b1}; // 동기 de-assert
    end
    assign rst_n_sync = shreg[STAGES-1];
endmodule
