`ifndef DIFFUSION_MATRIX_ADDKEY_V
`define DIFFUSION_MATRIX_ADDKEY_V


module MixColumns_AddKey
    (
        input [127:0] indata,
        input [127:0] key,
        output [127:0] outdata
    );
    localparam n = 128;
    localparam m = 4;
    genvar col, l;
    generate
        for(col=0; col< n/16; col=col+1) begin: gen_col 
            RotCol_AddKey rotcol_addkey(
            {indata[m*(col+1)-1:m*(col)], indata[m*(n/16+col+1)-1:m*(n/16+col)], indata[m*(n*2/16+col+1)-1:m*(n*2/16+col)], indata[m*(n*3/16+col+1)-1:m*(n*3/16+col)]}, 
                {key[m*(col+1)-1:m*(col)], key[m*(n/16+col+1)-1:m*(n/16+col)], key[m*(n*2/16+col+1)-1:m*(n*2/16+col)], key[m*(n*3/16+col+1)-1:m*(n*3/16+col)]},
                {outdata[m*(col+1)-1:m*(col)], outdata[m*(n/16+col+1)-1:m*(n/16+col)], outdata[m*(n*2/16+col+1)-1:m*(n*2/16+col)], outdata[m*(n*3/16+col+1)-1:m*(n*3/16+col)]}); 
        end
    endgenerate
endmodule



module RotCol_AddKey
    (
		inCols, key, outCols
    );
    localparam m = 4;
	input wire [m*4-1:0] inCols;
    input wire [m*4-1:0] key;
	output wire [m*4-1:0] outCols;
    genvar i;
    generate
        for(i=0; i<4; i=i+1)begin: gen_element
            wire [m*4-1:0] shiftedCol;
            if (i == 0) begin: branch1
                assign shiftedCol = inCols;
            end
            else begin: branch2
                assign shiftedCol = {inCols[(15+i*m)%16:0], inCols[15:i*m]};
            end
            assign outCols[m*(i+1)-1:i*m] = shiftedCol[2*m-1:m]^shiftedCol[3*m-1:2*m]^shiftedCol[4*m-1:3*m]^key[m*(i+1)-1:i*m];
        end
    endgenerate
endmodule

`endif