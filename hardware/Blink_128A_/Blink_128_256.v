`include "Blink_clk.v"
module Blink_128_256(clk, rst, enc, K0, P, T, C);
    input clk;
	input rst;
    localparam  round = 20; 
    localparam n = 128;
    localparam tweakLen = 128;
    localparam r1 = 4;
    localparam dim = 128;
    input enc;
	input [n*(round/2)-1:0] K0; 
    //input [509:0] K1; 
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