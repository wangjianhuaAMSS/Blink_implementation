`ifndef SHUFFLECELLS_V
`define SHUFFLECELLS_V
module ShuffleCells
    (
        input wire [63:0] indata,
        output wire [63:0] outdata
    );
    localparam n = 64;
	localparam m=4;
    genvar i;
    localparam [0:8*16-1] perm = {
    8'h00, 8'h05, 8'h0b, 8'h0a, 8'h01, 8'h06, 8'h04, 8'h0d,
    8'h02, 8'h0c, 8'h09, 8'h0f, 8'h03, 8'h07, 8'h0e, 8'h08
    };
    generate
        for(i=0; i< (n>>2); i=i+1)begin: gen_cell1
            assign outdata[(i+1)*m-1:i*m] = indata[(perm[i*8+:8]+1)*m-1:perm[i*8+:8]*m];
        end
    endgenerate
endmodule

module InvShuffleCells
    (
        input wire [63:0] indata,
        output wire [63:0] outdata
    );
    localparam n = 64;
	localparam m=4;
    genvar i;
    localparam [0:8*16-1] perm = {
    8'h00, 8'h05, 8'h0b, 8'h0a, 8'h01, 8'h06, 8'h04, 8'h0d,
    8'h02, 8'h0c, 8'h09, 8'h0f, 8'h03, 8'h07, 8'h0e, 8'h08
    };
    generate
        for(i=0; i< (n>>2); i=i+1)begin: gen_cell2
            assign outdata[(perm[i*8+:8]+1)*m-1:perm[i*8+:8]*m] = indata[(i+1)*m-1:i*m];
        end
    endgenerate
endmodule


`endif