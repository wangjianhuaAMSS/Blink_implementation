`ifndef SUBCELLS_V
`define SUBCELLS_V
`include "sbox.v"
module SubCells
    (
        input wire [127:0] indata,
        output wire [127:0] outdata
    );
	localparam m=32;
    genvar i;
    generate
        for(i=0; i<m; i=i+1)begin: gen_cell
            sbox  sbox_inst(indata[i*4+4-1:i*4], outdata[i*4+4-1:i*4]);
        end
    endgenerate
endmodule

`endif