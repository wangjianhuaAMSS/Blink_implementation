`include  "subcells.v"
`include  "diffusion_matrix_addkey.v"
`include "ShuffleCells.v"

module R
    (
        input wire [63:0] tk,
        input wire [63:0] indata,
        output wire [63:0] outdata
    );
        localparam blockLen = 64;
        wire [blockLen-1:0] first;

        SubCells sub_inst(.indata(indata), .outdata(first));
        wire [blockLen-1:0] middle;
        MixColumns_AddKey  mix_addkey_inst(.indata(first), .key(tk), .outdata(middle));
        ShuffleCells  shuffle_inst(.indata(middle), .outdata(outdata));
endmodule
