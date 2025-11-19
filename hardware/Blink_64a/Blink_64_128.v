`include "Blink_clk.v"
module Blink_64_128(clk, rst, enc, K0, P, T, C);
    input clk;
	input rst;
    localparam  round = 14; 
    localparam n = 64;
    localparam tweakLen = 64;
    localparam r1 = 3;
    localparam dim = 64;
    input enc;
	input [n*(round/2)-1:0] K0; 
    //input [(n*2-1)*2-1:0] K1; 
    input [n-1:0] P;
    input [tweakLen-1:0] T; 
    output [n-1:0] C;

    Blink_clk
    U_Blinkv2_clk
    (
        .clk(clk),
        .rst(rst),       
        .enc(enc),
        .K0(K0), 
        //.K1(K1), 
        .P(P), 
        .T(T), 
        .C(C)
    );
	
endmodule