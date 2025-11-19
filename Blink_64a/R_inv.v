`include  "subcells.v"
`include  "diffusion_matrix_addkey.v"
`include "ShuffleCells.v"

module R_inv
    (
        input wire [63:0] tk,
        input wire [63:0] indata,
        output wire [63:0] outdata
    );
    localparam blockLen = 64;
    wire [blockLen-1:0] first;
    wire [blockLen-1:0] middle;            
    InvShuffleCells  inv_shuffle_inst(.indata(indata), .outdata(first));
    MixColumns_AddKey  mix_addkey_inst(.indata(first), .key(tk), .outdata(middle));
    SubCells  sub_inst(.indata(middle), .outdata(outdata)); 
endmodule
