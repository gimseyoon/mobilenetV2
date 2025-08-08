module glbl_ctrl #(
    parameter IO_WIDTH = 18,
    parameter ROW = 14,
    parameter COLUMN = 14,
    parameter PIXEL = ROW * COLUMN,              // 14 * 14 = 196
    parameter W_WIDTH = 17,
    parameter ADDR_CHANNEL  = $clog2(384),        // 9 (for CHANNEL = 384)
    parameter ADDR_WMEM = $clog2(384 * 64),       // 15 (for 64*384 = 24576)
    parameter ADDR_W1_MEM = $clog2(384 * 9)       // 12 (for 9*384 = 3456)
)(
    input                               clk,
    input                               rst_n,
    input                        [2:0]  state,
    input                               pw_1_valid,
    input                               dw_valid,
    input                               save_valid,
    input           [ADDR_CHANNEL-1:0]  channel_num,
    input  signed [IO_WIDTH*PIXEL-1:0]  acc_out,
    input         [ADDR_CHANNEL -1 :0]  acc_cnt,
    
////////////////////////////////////////////////////////////////////////////
// bram_A
    output                              ena_0,
    output                              wea_0,
    output          [ADDR_CHANNEL-1:0]  addra_0,
    output reg signed [IO_WIDTH*PIXEL-1:0]  dina_0,
    output                              enb_0,
    output          [ADDR_CHANNEL-1:0]  addrb_0,

// bram_B
    output                              ena_1,
    output                              wea_1,
    output           [ADDR_CHANNEL-1:0] addra_1,
    output reg signed [IO_WIDTH*PIXEL-1:0]  dina_1,
    output                              enb_1,
    output           [ADDR_CHANNEL-1:0] addrb_1,

// bram_W0
    output                              ena_w0,
    output              [ADDR_WMEM-1:0] addra_w0,
// bram_W1
    output                              ena_w1,
    output              [ADDR_W1_MEM-1:0] addra_w1,
// bram_W2
    output                              ena_w2,
    output              [ADDR_WMEM-1:0] addra_w2,
////////////////////////////////////////////////////////////////////////////
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

    reg [14:0] glbl_cnt; // (0 ~ 32,767)


//////////////////////////////////////////////////////////////
// dina_0 / dina_1

    always @(*) begin
        if(!rst_n) begin
            dina_0 = 0; 
            dina_1 = 0;
        end
        else begin
            case (state)
                IDLE: begin
                    dina_0 = 0;
                    dina_1 = 0;
                end //IDLE
                
                PW_1: begin
                    if(pw_1_valid) begin
                        dina_1 = acc_out;
                        dina_0 = 0;
                    end  
                    else begin
                        dina_1 = 0;
                        dina_0 = 0;
                    end
                end //PW_1
                
                PW_1_BN_RELU: begin
                
                end
                
                DW: begin
                    if(dw_valid) begin
                        dina_1 = acc_out;
                        dina_0 = 0;
                    end  
                    else begin
                        dina_1 = 0;
                        dina_0 = 0;
                    end
                end //DW   

                // Other states: modify in future
                default: begin

                end
            endcase
        end
    end




//////////////////////////////////////////////////////////////
// cnt : 0 -> 24579 -> 0 -> 24579 -> ...

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        glbl_cnt <= 0;
    end else begin
        case (state)
            PW_1: begin
                if (glbl_cnt >= 15'd24576 + 15'd7)
                    glbl_cnt <= 0;
                else
                    glbl_cnt <= glbl_cnt + 1;
            end
            
            PW_1_BN_RELU: begin
                
            end
            
            DW: begin
                if (glbl_cnt >= 15'd3256 + 15'd3)
                    glbl_cnt <= 0;
                else
                    glbl_cnt <= glbl_cnt + 1;
            end
            default: begin
                glbl_cnt <= 0;
            end
        endcase
    end
end

//////////////////////////////////////////////////////////////
// instantiation addr_counter

    addr_counter addr_counter_0 (
        .clk                (clk),
        .rst_n              (rst_n),
        .state              (state),
        .glbl_cnt           (glbl_cnt),
        .acc_cnt            (acc_cnt),
        .save_valid         (save_valid),
        .pw_1_valid         (pw_1_valid),
        .channel_num        (channel_num),
        .enb_0              (enb_0),
        .enb_1              (enb_1),
        
        .addra_0           (addra_0),
        .addrb_0           (addrb_0),
        .addra_1           (addra_1),
        .addrb_1           (addrb_1),
        .addra_w0          (addra_w0),
        .addra_w1          (addra_w1),
        .addra_w2          (addra_w2),
        .addra_bias_0      (addra_bias_0),
        .addra_mean_0      (addra_mean_0),
        .addra_std_0       (addra_std_0),
        .addra_weight_0    (addra_weight_0)
    );

//////////////////////////////////////////////////////////////
// instantiation enable_counter
    
    enable_counter enable_counter_0 (
        .clk              (clk),
        .rst_n            (rst_n),
        .state            (state),
        .glbl_cnt         (glbl_cnt),
        .acc_cnt          (acc_cnt),
        .pw_1_valid       (pw_1_valid),
        
        .ena_0           (ena_0),
        .wea_0           (wea_0),
        .enb_0           (enb_0),
        .ena_1           (ena_1),
        .wea_1           (wea_1),
        .enb_1           (enb_1),
        .ena_w0          (ena_w0),
        .ena_w1          (ena_w1),
        .ena_w2          (ena_w2),
        .ena_bias_0      (ena_bias_0),
        .ena_mean_0      (ena_mean_0),
        .ena_std_0       (ena_std_0),
        .ena_weight_0    (ena_weight_0)
    );

endmodule