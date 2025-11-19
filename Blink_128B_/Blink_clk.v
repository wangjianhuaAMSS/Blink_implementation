`include "Blink_top.v"
module Blink_clk
	(clk, rst, enc, K0, P, T, C);
    input clk;
	input rst;
    localparam n = 128;
    localparam tweakLen = 256;
    localparam round = 20;
    localparam dim = 128;

    input enc;
    input [n*(round/2)-1:0] K0;
    //input [(128*3-1)*2-1:0] K1;
    input [n-1:0] P;
    input [tweakLen-1:0] T;
    output [n-1:0] C;
    reg [n-1:0] C;

    wire [n-1:0] C_wire;

    reg [n*(round/2)-1:0] K0_reg;
    //reg [(128*3-1)*2-1:0] K1_reg; 
    reg [n-1:0] P_reg;
    reg [tweakLen-1:0] T_reg; 

    Blink_top
    Blinkv2_top_inst
    (   
        .enc(enc),
        .K0(K0_reg), 
        //.K1(K1_reg), 
        .P(P_reg), 
        .T(T_reg), 
        .C(C_wire)
    );

	always@(posedge clk or negedge rst) begin
        if (rst == 0) begin: branch1
            K0_reg <= 0;
            //K1_reg <= 0; 
            P_reg <= 0;
            T_reg <= 0; 
            C <= 0;
        end
        else begin: branch2
            K0_reg <= K0;
            //K1_reg <= K1; 
            P_reg <= P;
            T_reg <= T; 
            C <= C_wire;
        end

    end
	
endmodule