`timescale 1ns / 1ps

module accumulator#(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN, // // 14 * 14 = 196 PIXEL
    parameter W_WIDTH = 17
)(
    input                                      clk,
    input                                      rst_n,
    input                              [2:0]   state,
    input  signed [IO_WIDTH * PIXEL - 1 : 0]   mul_out,  // [3920-1 : 0], 3920 bit
    output reg                                 save_valid,
    output reg                         [8:0]   channel_num,
    output reg                                 pw_1_done,
    output signed [IO_WIDTH * PIXEL - 1 : 0]   acc_out    // [3920-1 : 0], 3920 bit
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

    reg signed [IO_WIDTH-1:0] acc_out_reg [0:PIXEL-1];  // [20-1: 0]] acc_out_reg [196-1 :0]*/
    reg [8:0] cnt; // Cover 0 ~ 383
    reg [7:0] state_delay;
    
    wire conv_en;
    wire signed [IO_WIDTH-1:0] mul_out_array [0:PIXEL-1]; // Convert [3920-1 :0] mul_out -> [20-1 :0] mul_out_array [196-1 :0]

///////////////////////////////////////////////////////////////////////
// assign conv_en = state_delay[7]

    assign conv_en = state_delay[7];



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
            if(state == PW_1) begin
                state_delay[0] <= 1;
            end
            else begin
                state_delay <=0 ;
            end
        end //else
    end //always

///////////////////////////////////////////////////////////////////////
// assign mul_out_array[i] = mul_out[20(i+1) -1 :20*i]

    genvar j;
    generate
      for (j = 0; j < PIXEL; j = j + 1) begin : UNPACK_MUL_OUT
        assign mul_out_array[j] = mul_out[IO_WIDTH*(j+1)-1 : IO_WIDTH*j];   // Convert [3920-1 :0] mul_out -> [20-1 :0] mul_out_array [196-1 :0]
      end
    endgenerate
    
    
    
///////////////////////////////////////////////////////////////////////
// assign acc_out[3920-1 :0] = acc_out_reg [20-1 :0][196-1 :0]

    genvar i;
    generate
        for (i = 0; i < PIXEL; i = i + 1) begin : PACK_OUTPUT
            assign acc_out[IO_WIDTH*(i+1)-1 : IO_WIDTH*i] = (save_valid) ? acc_out_reg[i] : 0;
        end
    endgenerate
    
    


///////////////////////////////////////////////////////////////////////
// pointwise's output, depthwise's output, save_valid signal

    integer k;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (k = 0; k < PIXEL; k = k + 1) begin
                acc_out_reg[k] <= 0;
            end
            cnt <= 0;
            save_valid <= 0;
            channel_num <= 0;
            pw_1_done <= 0;
        end else begin
            case (state)
                IDLE: begin
                    for (k = 0; k < PIXEL; k = k + 1) begin
                        acc_out_reg[k] <= 0;
                    end
                    cnt <= 0;
                    save_valid <= 0;
                    channel_num <= 0;
                end
    
                PW_1: begin
                
                    for (k = 0; k < PIXEL; k = k + 1) begin
                        if (cnt >= 63)  begin acc_out_reg[k] <= mul_out_array[k];                  end // when start a new channel
                        else            begin acc_out_reg[k] <= acc_out_reg[k] + mul_out_array[k]; end //accumulate 
                    end
                    
                    if (conv_en) begin
                        cnt <= (cnt == 63) ? 0 : cnt + 1;
                        save_valid <= (cnt == 62);
                        pw_1_done <= (cnt==62) && (channel_num == 383);
                        
                        if (save_valid) begin channel_num <= (channel_num == 383) ? 0 : channel_num + 1; end
                        
                    end
                    else begin
                        // conv_en ¨¬??¡Æ¨ù¨¬ ¨ö? ??¡¾??¡©
                        cnt <= 0;
                        save_valid <= 0;
                        channel_num <= 0;
                        pw_1_done <= 0;
                    end
                end // PW_1

                
                PW_1_BN_RELU:
                    begin
                    
                    end
    
                DW: begin

                end
                
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
                    cnt <= 0;
                    save_valid <= 0;
                    channel_num <= 0;
                    pw_1_done <= 0;
                end
            endcase
        end
    end

endmodule