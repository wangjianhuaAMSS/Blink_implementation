`ifndef SHUFFLECELLS_V
`define SHUFFLECELLS_V
module ShuffleCells
    (
        input wire [127:0] indata,
        output wire [127:0] outdata
    );
    localparam n = 128;
	localparam m=4;
    genvar i;
    localparam [0:8*32-1] perm = {
    8'h05, 8'h0c, 8'h04, 8'h01, 8'h11, 8'h09, 8'h0a, 8'h10, 
    8'h1c, 8'h0e, 8'h15, 8'h16, 8'h0b, 8'h1b, 8'h08, 8'h0d, 
    8'h02, 8'h19, 8'h12, 8'h03, 8'h1e, 8'h06, 8'h13, 8'h14, 
    8'h00, 8'h17, 8'h18, 8'h1f, 8'h07, 8'h0f, 8'h1d, 8'h1a
    };

    generate
        for(i=0; i< (n>>2); i=i+1)begin: gen_cell1
            assign outdata[(i+1)*m-1:i*m] = indata[(perm[i*8+:8]+1)*m-1:perm[i*8+:8]*m];
        end
    endgenerate
endmodule

module InvShuffleCells
    (
        input wire [127:0] indata,
        output wire [127:0] outdata
    );
    localparam n = 128;
	localparam m=4;
    genvar i;

    localparam [0:8*32-1] perm = {
    8'h05, 8'h0c, 8'h04, 8'h01, 8'h11, 8'h09, 8'h0a, 8'h10, 
    8'h1c, 8'h0e, 8'h15, 8'h16, 8'h0b, 8'h1b, 8'h08, 8'h0d, 
    8'h02, 8'h19, 8'h12, 8'h03, 8'h1e, 8'h06, 8'h13, 8'h14, 
    8'h00, 8'h17, 8'h18, 8'h1f, 8'h07, 8'h0f, 8'h1d, 8'h1a
    };
    generate
        for(i=0; i< (n>>2); i=i+1)begin: gen_cell2
            assign outdata[(perm[i*8+:8]+1)*m-1:perm[i*8+:8]*m] = indata[(i+1)*m-1:i*m];
        end
    endgenerate
endmodule


`endif