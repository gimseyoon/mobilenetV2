module glbl_ctrl #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),         // 8 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64)       // 15 (for 64*384 = 24576)
)(
    input                               clk,
    input                               rst_n,
    input                        [2:0]  state,
    input                               save_valid,
    input           [ADDR_CHANNEL-1:0]  channel_num,
    input                               bram_select,
    input  signed [IO_WIDTH*PIXEL-1:0]  acc_out,
    input         [ADDR_CHANNEL -1 :0]  acc_cnt,
    
// bram_0
    output                              ena_0,
    output                              wea_0,
    output          [ADDR_CHANNEL-1:0]  addra_0,
    output signed [IO_WIDTH*PIXEL-1:0]  dina_0,
    output                              enb_0,
    output          [ADDR_CHANNEL-1:0]  addrb_0,

// bram_1
    output                              ena_1,
    output                              wea_1,
    output           [ADDR_CHANNEL-1:0] addra_1,
    output signed [IO_WIDTH*PIXEL-1:0]  dina_1,
    output                              enb_1,
    output           [ADDR_CHANNEL-1:0] addrb_1,

// bram_w
    output                              ena_w0,
    output              [ADDR_WMEM-1:0] addra_w0,
    
// BRAM bias
    output                              ena_bias_0,    
    output         [ADDR_CHANNEL-1 : 0] addra_bias_0,
// BRAM mean
    output                              ena_mean_0,    
    output         [ADDR_CHANNEL-1 : 0] addra_mean_0,  
// BRAM std
    output                              ena_std_0,     
    output         [ADDR_CHANNEL-1 : 0] addra_std_0,   
// BRAM weight    
    output                              ena_weight_0,  
    output         [ADDR_CHANNEL-1 : 0] addra_weight_0
);

//////////////////////////////////////////////////////////////

    localparam IDLE         = 3'b000,
               PW_1         = 3'b001,
               PW_1_BN_RELU = 3'b010,
               DW           = 3'b011,
               DW_BN_RELU   = 3'b100,
               PW_2         = 3'b101,
               PW_2_BN      = 3'b110,
               SK           = 3'b111;

//////////////////////////////////////////////////////////////

    reg [14:0] cnt;

    // input bram
    wire in_ena;
    wire in_wea;
    wire [ADDR_CHANNEL-1 : 0] in_addra;
    wire in_enb;
    wire [ADDR_CHANNEL-1 : 0] in_addrb;
    // output bram
    wire out_ena;
    wire out_wea;
    wire [ADDR_CHANNEL-1 : 0] out_addra;
    wire out_enb;
    wire [ADDR_CHANNEL-1 : 0] out_addrb;
    //weight bram
    wire weight_ena;
    wire [ADDR_WMEM-1 : 0] weight_addra;
    
//////////////////////////////////////////////////////////////

// [bram_0] : bram_select == 1 -> INPUT, bram_select == 0 -> OUTPUT 
    assign ena_0     = (bram_select == 1) ? 0          : save_valid;
    assign wea_0     = (bram_select == 1) ? 0          : save_valid;
    assign addra_0   = (bram_select == 1) ? 0          : channel_num;
    assign dina_0    = (bram_select == 1) ? 0          : acc_out;
    assign enb_0     = (bram_select == 1) ? in_enb     : 0;
    assign addrb_0   = (bram_select == 1) ? in_addrb   : 0;
    
// [bram_1] : bram_select == 0 -> INPUT, bram_select == 1 -> OUTPUT 
    assign ena_1     = (bram_select == 0) ? 0          : save_valid;
    assign wea_1     = (bram_select == 0) ? 0          : save_valid;
    assign addra_1   = (bram_select == 0) ? 0          : channel_num;
    assign dina_1    = (bram_select == 0) ? 0          : acc_out;
    assign enb_1     = (bram_select == 0) ? in_enb     : 0;
    assign addrb_1   = (bram_select == 0) ? in_addrb   : 0;

// bram_w0
    assign ena_w0    = weight_ena;
    assign addra_w0  = weight_addra;


//////////////////////////////////////////////////////////////
// cnt : 0 -> 24579 -> 0 -> 24579 -> ...

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            cnt <= 0;
        end else begin
            if (state == PW_1) begin
                if (cnt >= 15'd24576 + 15'd3)
                    cnt <= 0;
                else
                    cnt <= cnt + 1;
            end else begin
                cnt <= 0;
            end
        end
    end

//////////////////////////////////////////////////////////////
// instantiation addr_counter

    addr_counter addr_counter_0 (
        .clk(clk),
        .rst_n(rst_n),
        .state(state),
        .cnt(cnt),
        .acc_cnt(acc_cnt),
        .save_valid(save_valid),
        .channel_num(channel_num),
        .in_ena(in_ena),
        .in_enb(in_enb),
        .weight_ena(weight_ena),
        .in_addra(in_addra),
        .in_addrb(in_addrb),
        .weight_addra(weight_addra)
    );

//////////////////////////////////////////////////////////////
// instantiation enable_counter
    
    enable_counter enable_counter_0 (
        .clk(clk),
        .rst_n(rst_n),
        .state(state),
        .cnt(cnt),
        .acc_cnt(acc_cnt),
        .in_ena(in_ena),
        .in_wea(in_wea),
        .in_enb(in_enb),
        .weight_ena(weight_ena)
    );

endmodule