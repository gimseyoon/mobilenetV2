`timescale 1ns / 1ps

module accumulator #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),        // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64),       // 15 (for 64*384 = 24576)
    parameter ADDR_W1_MEM = $clog2(384 * 9)       // 12 (for 9*384 = 3456)
)(
    input                                      clk,
    input                                      rst_n,
    input                              [2:0]   state,
    input  signed  [IO_WIDTH * PIXEL - 1 :0]   mul_out,  // [3920-1 : 0], 3920 bit
    output reg                                 bn_en,
    output reg         [ADDR_CHANNEL -1 :0]    acc_cnt,
    output reg                                 pw_1_valid,
    output reg                                 pw_1_done,
    output reg                                 dw_valid,
    output reg                                 dw_done,
    output reg                                 pw_2_valid,
    output reg                                 pw_2_done,
    output reg          [ADDR_CHANNEL-1 :0]    channel_num,
    output signed [IO_WIDTH * PIXEL - 1 : 0]    acc_out    // [3920-1 : 0], 3920 bit
);

////////////////////////////////////////////////////////////

    localparam IDLE         = 3'b000,
               PW_1         = 3'b001,
               PW_1_BN_RELU = 3'b010,
               DW           = 3'b011,
               DW_BN_RELU   = 3'b100,
               PW_2         = 3'b101,
               PW_2_BN      = 3'b110,
               SK           = 3'b111;

////////////////////////////////////////////////////////////

    reg [3:0] bn_en_cnt;
    reg signed [IO_WIDTH-1:0] acc_out_reg [0:PIXEL-1];  // [20-1: 0]] acc_out_reg [196-1 :0]*/
    reg [7:0] state_delay;
    reg signed [IO_WIDTH-1:0] temp [0:PIXEL-1];
    wire pw_1_en;
    wire dw_en;
    wire signed [IO_WIDTH-1:0] mul_out_reg [0:PIXEL-1]; // Convert [3920-1 :0] mul_out -> [20-1 :0] mul_out_reg [196-1 :0]

///////////////////////////////////////////////////////////////////////
// assign pw_1_en = state_delay[7]

    assign pw_1_en = (state==PW_1) ? state_delay[7] : 0;
    assign dw_en = (state==DW) ? state_delay[5] : 0;


///////////////////////////////////////////////////////////////////////
// state_delay (to wait bram read)

    always@(posedge clk or negedge rst_n) begin
        if(!rst_n) begin
           state_delay <= 0; 
        end
        else begin
            state_delay[7] <= state_delay[6];
            state_delay[6] <= state_delay[5];
            state_delay[5] <= state_delay[4];
            state_delay[4] <= state_delay[3];
            state_delay[3] <= state_delay[2];
            state_delay[2] <= state_delay[1];
            state_delay[1] <= state_delay[0];
            if(pw_1_done || dw_done || pw_2_done) begin
                state_delay <= 0;
            end
            else if( (state == PW_1) || (state == DW) || (state == PW_2)) begin
                state_delay[0] <= 1;
            end
            else begin
                state_delay <=0 ;
            end
        end //else
    end //always

///////////////////////////////////////////////////////////////////////
// assign mul_out_reg[i] = mul_out[20(i+1) -1 :20*i]

    genvar j;
    generate
      for (j = 0; j < PIXEL; j = j + 1) begin : UNPACK_MUL_OUT
        assign mul_out_reg[j] = mul_out[IO_WIDTH*(j+1)-1 : IO_WIDTH*j];   // Convert [3920-1 :0] mul_out -> [20-1 :0] mul_out_reg [196-1 :0]
      end
    endgenerate
    
    
    
///////////////////////////////////////////////////////////////////////
// assign acc_out[3920-1 :0] = acc_out_reg [20-1 :0][196-1 :0]

    genvar i;
    generate
        for (i = 0; i < PIXEL; i = i + 1) begin : PACK_OUTPUT
            assign acc_out[IO_WIDTH*(i+1)-1 : IO_WIDTH*i] = (pw_1_valid || dw_valid || pw_2_valid) ? acc_out_reg[i] : 0;
        end
    endgenerate
    
    


///////////////////////////////////////////////////////////////////////
// pointwise's output, depthwise's output, bn_relu_valid signal

    integer k, x, y, p;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (k = 0; k < PIXEL; k = k + 1) begin
                acc_out_reg[k] <= 0;
            end
            bn_en_cnt <= 0;
            acc_cnt <= 0;
            channel_num <= 0;
            bn_en <= 0;
            pw_1_valid <= 0;
            pw_1_done <= 0;
            dw_valid <= 0;
            dw_done <= 0;
            pw_2_valid <= 0;
            pw_2_done <= 0;
        end 
        else begin
            case (state)
                IDLE: begin
                    for (k = 0; k < PIXEL; k = k + 1) begin
                        acc_out_reg[k] <= 0;
                    end
                    bn_en_cnt <= 0;
                    acc_cnt <= 0;
                    channel_num <= 0;
                    bn_en <= 0;
                    pw_1_valid <= 0;
                    pw_1_done <= 0;
                    dw_valid <= 0;
                    dw_done <= 0;
                    pw_2_valid <= 0;
                    pw_2_done <= 0;
                end
    
                PW_1: begin
                
                    for (k = 0; k < PIXEL; k = k + 1) begin
                        if (acc_cnt >= 63)  begin acc_out_reg[k] <= mul_out_reg[k];                  end // when start a new channel
                        else                begin acc_out_reg[k] <= acc_out_reg[k] + mul_out_reg[k]; end //accumulate 
                    end
                    
                    if (pw_1_en) begin
                        acc_cnt <= (acc_cnt == 63) ? 0 : acc_cnt + 1;
                        pw_1_valid <= (acc_cnt == 62);
                        pw_1_done <= (acc_cnt==62) && (channel_num == 383);
                        
                        if (pw_1_valid) begin channel_num <= (channel_num == 383) ? 0 : channel_num + 1; end
                    end
                    else begin
                        acc_cnt <= 0;
                        pw_1_valid <= 0;
                        pw_1_done <= 0;
                        channel_num <= 0;
                    end
                    
                    if(pw_1_done) begin
                        for (k = 0; k < PIXEL; k = k + 1) begin
                            acc_out_reg[k] <= 0;
                        end
                    end
                    
                    if(pw_1_valid) begin
                        bn_en <= 1;
                    end
                    else begin
                        if(bn_en_cnt == 13) bn_en <= 0;
                    end
                    
                    if(bn_en) begin
                        if(bn_en_cnt == 13) bn_en_cnt <= 0;
                        else                bn_en_cnt <= bn_en_cnt + 1;
                    end
                    else begin
                        bn_en_cnt <= 0;
                    end
                    
                    
                    
                end // PW_1

                
                PW_1_BN_RELU:
                    begin
                    
                    end
    
                DW: begin
                
                        bn_en <= dw_valid && (channel_num != 64);   
                        
                        case (acc_cnt)
                                // case 1: ?????? ??? 13x13
                                0: begin
                                    for (y = 0; y < 13; y = y + 1) begin
                                        for (x = 0; x < 13; x = x + 1) begin
                                            acc_out_reg[(y + 1) * 14 + (x + 1)] <= mul_out_reg[y * 14 + x];
                                        end
                                    end
                                end
                    
                                // case 2: ????? 13x14
                                1: begin
                                    for (y = 0; y < 13; y = y + 1) begin
                                        for (x = 0; x < 14; x = x + 1) begin
                                            acc_out_reg[(y + 1) * 14 + x] <= acc_out_reg[(y + 1) * 14 + x] + mul_out_reg[y * 14 + x];
                                        end
                                    end
                                end
                    
                                // case 3: ???? ??? 13x13
                                2: begin
                                    for (y = 0; y < 13; y = y + 1) begin
                                        for (x = 1; x < 14; x = x + 1) begin
                                            acc_out_reg[(y + 1) * 14 + (x - 1)] <= acc_out_reg[(y + 1) * 14 + (x - 1)] + mul_out_reg[y * 14 + x];
                                        end
                                    end
                                end
                    
                                // case 4: ?????? 14x13
                                3: begin
                                    for (y = 0; y < 14; y = y + 1) begin
                                        for (x = 0; x < 13; x = x + 1) begin
                                            acc_out_reg[y * 14 + (x + 1)] <= acc_out_reg[y * 14 + (x + 1)] + mul_out_reg[y * 14 + x];
                                        end
                                    end
                                end
                    
                                // case 5: ??? 14x14 ????
                                4: begin
                                    for (p = 0; p < 196; p = p + 1) begin
                                        acc_out_reg[p] <= acc_out_reg[p] + mul_out_reg[p];
                                    end
                                end
                    
                                // case 6: ???? 14x13
                                5: begin
                                    for (y = 0; y < 14; y = y + 1) begin
                                        for (x = 1; x < 14; x = x + 1) begin
                                            acc_out_reg[y * 14 + (x - 1)] <= acc_out_reg[y * 14 + (x - 1)] + mul_out_reg[y * 14 + x];
                                        end
                                    end
                                end
                    
                                // case 7: ?????? ?? 13x13
                                6: begin
                                    for (y = 1; y < 14; y = y + 1) begin
                                        for (x = 0; x < 13; x = x + 1) begin
                                            acc_out_reg[(y - 1) * 14 + (x + 1)] <= acc_out_reg[(y - 1) * 14 + (x + 1)] + mul_out_reg[y * 14 + x];
                                        end
                                    end
                                end
                    
                                // case 8: ???? 13x14
                                7: begin
                                    for (y = 1; y < 14; y = y + 1) begin
                                        for (x = 0; x < 14; x = x + 1) begin
                                            acc_out_reg[(y - 1) * 14 + x] <= acc_out_reg[(y - 1) * 14 + x] + mul_out_reg[y * 14 + x];
                                        end
                                    end
                                end
                    
                                // case 9: ???? ?? 13x13
                                8: begin
                                    for (y = 1; y < 14; y = y + 1) begin
                                        for (x = 1; x < 14; x = x + 1) begin
                                            acc_out_reg[(y - 1) * 14 + (x - 1)] <= acc_out_reg[(y - 1) * 14 + (x - 1)] + mul_out_reg[y * 14 + x];
                                        end
                                    end
                                end
                                
                                9: begin
                                
                                    // 1. ??? ???? ????
                                    for (k = 0; k < PIXEL; k = k + 1) begin
                                        temp[k] = 0;
                                    end
                                
                                    // 2. ?????? ??? 13??13 ????
                                    for (y = 0; y < 13; y = y + 1) begin
                                        for (x = 0; x < 13; x = x + 1) begin
                                            temp[(y + 1) * 14 + (x + 1)] = mul_out_reg[y * 14 + x];
                                        end
                                    end
                                
                                    // 3. acc_out_reg?? ???
                                    for (k = 0; k < PIXEL; k = k + 1) begin
                                        acc_out_reg[k] <= temp[k];
                                    end
                                end
                                                                
                                
                                // default (optional): ??? ????? ???? ????
                                default: begin
                                    // no operation
                                end
                            endcase
                    
                    if (dw_en) begin
                        acc_cnt <= (acc_cnt == 9) ? 1 : acc_cnt + 1;
                        dw_valid <= (acc_cnt == 8);
                        dw_done <= (acc_cnt==8) && (channel_num == 383);
                        
                        if (dw_valid) begin channel_num <= (channel_num == 383) ? 0 : channel_num + 1; end
                    end
                    else begin
                        acc_cnt <= 511;
                        dw_valid <= 0;
                        dw_done <= 0;
                        channel_num <= 0;
                    end
                    
                    if(dw_done) begin
                        for (k = 0; k < PIXEL; k = k + 1) begin
                            acc_out_reg[k] <= 0;
                        end                    
                    end
                    
                end //DW
                
                DW_BN_RELU:
                    begin
                    
                    end
                    
                PW_2: begin

                end
                
                PW_2_BN:
                    begin
                    
                    end
                    
                SK:
                    begin
                    
                    end
                
                default: begin
                    for (k = 0; k < PIXEL; k = k + 1) begin
                        acc_out_reg[k] <= 0;
                    end
                    acc_cnt <= 0;
                    channel_num <= 0;
                    pw_1_valid <= 0;
                    pw_1_done <= 0;
                    dw_valid <= 0;
                    dw_done <= 0;
                    pw_2_valid <= 0;
                    pw_2_done <= 0;
                end
            endcase
        end
    end

endmodule