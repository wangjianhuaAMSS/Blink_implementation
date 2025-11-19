module generate_tweakey (tweak, key, tweakey);
    localparam blockLen = 64;
    localparam dim = 64;
    localparam tweakLen = 128;
    localparam TblockNum = blockLen / dim;
    localparam hashInputLen = tweakLen / TblockNum;
    localparam hashKLen = dim + hashInputLen - 1;
    input [tweakLen-1: 0] tweak;
    input [hashKLen*TblockNum*2-1: 0] key;
    output [blockLen * 2-1:0] tweakey;
    genvar i;

    generate
      for (i = 0; i <= TblockNum-1; i = i+1) begin: gen_tk02
        hash  H1 (.f_in(tweak[hashInputLen*(i+1)-1:hashInputLen*i]), .key(key[hashKLen*(i+1)-1:hashKLen*i]), .hash(tweakey[dim*(i+1)-1: dim *i]));
        hash  H2 (.f_in(tweak[hashInputLen*(i+1)-1:hashInputLen*i]), .key(key[hashKLen*(i+1)+hashKLen*TblockNum -1:hashKLen*i+hashKLen*TblockNum]), .hash(tweakey[dim*(i+1)-1+blockLen: dim *i+blockLen]));
      end
    endgenerate
endmodule

module hash (f_in, key, hash); 
  localparam inputLen = 128;
  localparam dim = 64; 
  input [inputLen-1:0] f_in;
  input [inputLen+dim-2: 0] key;
  output  [dim - 1: 0] hash;
  genvar i, k;

  generate
    for (k = 0; k < dim; k = k+1) begin: gen_k
      wire [dim-1:0] ThreeSum [2:0];
      innerproduct  IP (.f_in0(key[inputLen + dim - 2 - k: dim - 1 - k]), .f_in1(f_in), .f_out(hash[k]));
    end
  endgenerate
endmodule



module innerproduct  (f_in0, f_in1, f_out);
  localparam len = 128;
  input [len-1: 0] f_in0;
  input [len-1: 0] f_in1;
  output f_out;
  wire [len-1:0] tmp;
  genvar i;
  assign tmp[0] = f_in0[0] & f_in1[0];
  generate 
    for (i = 1; i < len; i =i+1) begin: gen_i
      assign tmp[i] = (f_in0[i] & f_in1[i]) ^ tmp[i - 1];
    end
  endgenerate
  assign f_out = tmp[len-1];

endmodule