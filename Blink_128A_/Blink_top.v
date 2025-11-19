`include "R.v"
`include "R_inv.v"
`include "generate_tweakey.v"
`include "diffusion_matrix.v"

module Blink_top
    (        
        input wire enc,
        input wire[128*10-1:0] K0, 
        //input wire[(128*2-1)*2-1:0] K1, 
        input wire[128-1:0] P, 
        input wire[128-1:0] T, 
        output wire[128-1:0] C
    );
    localparam [128*8-1:0] rc1 = { 
    128'hb4cc5c341141e8cea15486af7c72e993,
    128'h5748986263e8144055ca396a2aab10b6,
    128'h78af2fda55605c60e65525f3aa55ab94,
    128'h6c9e0e8bb01e8a3ed71577c1bd314b27,
    128'hca417918b8db38ef8e79dcb0603a180e,
    128'h9c30d5392af26013c5d1b023286085f0,
    128'h718bcd5882154aee7b54a41dc25a59b5,
    128'ha458fea3f4933d7e0d95748f728eb658};
    localparam [128*8-1:0] rc0 = {
    128'h0801f2e2858efc16636920d871574e69,
    128'hba7c9045f12c7f9924a19947b3916cf7,
    128'h2ffd72dbd01adfb7b8e1afed6a267e96,
    128'h9216d5d98979fb1bd1310ba698dfb5ac,
    128'hc0ac29b7c97c50dd3f84d5b5b5470917,
    128'h452821e638d01377be5466cf34e90c6c,
    128'ha4093822299f31d0082efa98ec4e6c89,
    128'h243f6a8885a308d313198a2e03707344
    };
    localparam  round = 20; 
    localparam blockLen = 128;
    localparam tweakLen = 128;
    localparam r1 = 4;

    genvar i;

    //tweakey generation
    wire [blockLen*(round/2-2)-1:0] mixrk;
    wire [blockLen*(round/2-2)-1:0] mixrc0;
    wire [blockLen*(round/2-2)-1:0] mixrc1;
    wire [blockLen*2-1:0] mixhash;
    wire [blockLen*(round-1)-1:0] tweakey;
    wire [2*blockLen-1: 0]  hash;
    wire [blockLen-1: 0]  wk[1:0];

    assign wk[0] = K0[blockLen-1:0];
    assign wk[1] = K0[2*blockLen-1:blockLen];

    wire[(128*2-1)*2-1:0] K1;
    generate
        for (i = 0; i <= (128*2-1)*2-1; i =i+1 ) begin:gen_k1
            assign K1[i] = K0[11*i%(128*10)];
        end
    endgenerate

    generate_tweakey  GTK (.tweak(T), .key(K1), .tweakey(hash));
    generate
        for (i = 1; i < round/2; i = i+1) begin: gen_tweakeyTmp
            if (i == r1) begin: branch11
                MixColumns  mixh1 (.indata(hash[blockLen-1: 0]), .outdata(mixhash[blockLen-1:0]));
                MixColumns  mixh2 (.indata(hash[2*blockLen-1: blockLen]), .outdata(mixhash[2*blockLen-1:blockLen]));
            end
            else begin: branch12
                MixColumns  mixRK (.indata(K0[(i+2-(i>r1))*blockLen-1:(i+1-(i>r1))*blockLen]), .outdata(mixrk[(i-(i>r1))*blockLen-1:(i-(i>r1)-1)*blockLen]));
                MixColumns  mixc0 (.indata(rc0[(i-(i>r1)-1)*128+blockLen-1: (i-(i>r1)-1)*128]), .outdata(mixrc0[(i-(i>r1))*blockLen-1:(i-(i>r1)-1)*blockLen]));
                MixColumns  mixc1 (.indata(rc1[(i-(i>r1)-1)*128+blockLen-1: (i-(i>r1)-1)*128]), .outdata(mixrc1[(i-(i>r1))*blockLen-1:(i-(i>r1)-1)*blockLen]));
            end   
        end
        for (i = 1; i <= round/2; i = i+1) begin: gen_tweakey
            if (i == r1) begin: branch21
                assign tweakey[i*blockLen-1:(i-1)*blockLen] = enc? hash[blockLen-1: 0]: hash[2*blockLen-1: blockLen];
                assign tweakey[(round-i)*blockLen-1:(round-i-1)*blockLen] = enc? mixhash[2*blockLen-1: blockLen]: mixhash[blockLen-1: 0];
            end 
            else if (i == round/2) begin: branch22
                assign tweakey[i*blockLen-1:(i-1)*blockLen] = enc? hash[blockLen-1: 0] ^ hash[2*blockLen-1: blockLen]: mixhash[blockLen-1: 0] ^ mixhash[2*blockLen-1: blockLen];
            end
            else begin: branch23
                assign tweakey[i*blockLen-1:(i-1)*blockLen] = enc? K0[(i+2-(i>r1))*blockLen-1:(i+1-(i>r1))*blockLen] ^ rc0[(i-(i>r1)-1)*128+blockLen-1: (i-(i>r1)-1)*128]: K0[(round/2-i+1+(i>r1))*blockLen-1:(round/2-i+(i>r1))*blockLen] ^ rc1[(round/2-2-i+(i>r1))*128+blockLen-1: (round/2-i+(i>r1)-2)*128];

                assign tweakey[(round-i)*blockLen-1:(round-i-1)*blockLen] = enc? mixrk[(round/2-1-(i-(i>r1)))*blockLen-1: (round/2-2-(i-(i>r1)))*blockLen] ^ mixrc1[(round/2-1-(i-(i>r1)))*blockLen-1: (round/2-1-(i-(i>r1))-1)*blockLen]: mixrk[(i-(i>r1))*blockLen-1: (i-(i>r1)-1)*blockLen] ^ mixrc0[(i-(i>r1))*blockLen-1: (i-(i>r1)-1)*blockLen];
            end

            
        end
    endgenerate 

    //data
    wire [blockLen-1:0] tmpd[round-1:0];
    
    assign tmpd[0] = enc? P ^ wk[0]:P ^ wk[1];

    generate
        for (i = 1; i <= round/2-1; i = i+1) begin: gen_rf
            R 
            rf
            (
                .tk(tweakey[(i)*blockLen-1:(i-1)*blockLen]),
                .indata(tmpd[i-1]),
                .outdata(tmpd[i])
            );
            R_inv 
            rf_inv
            (
                .tk(tweakey[(i+round/2)*blockLen-1:(i+round/2-1)*blockLen]),
                .indata(tmpd[round/2+i-1]),
                .outdata(tmpd[round/2+i])
            );
        end
        wire [blockLen-1:0] first;
        SubCells  sub_inst1(.indata(tmpd[round/2-1]), .outdata(first));
        wire [blockLen-1:0] middle;
        MixColumns_AddKey  mix_addkey_inst(.indata(first), .key(tweakey[(round/2)*blockLen-1:(round/2-1)*blockLen]), .outdata(middle));
        SubCells sub_inst2(.indata(middle), .outdata(tmpd[round/2]));
    endgenerate
    assign C = enc? tmpd[round-1] ^ wk[1]: tmpd[round-1] ^ wk[0];
endmodule
