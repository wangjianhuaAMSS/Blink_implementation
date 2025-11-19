#include <iostream>
#include <iomanip>
#include <stdio.h>
using namespace std;

constexpr int STATE_BYTES = 16;
constexpr int STATE_NIBBLES = STATE_BYTES * 2;
constexpr int TWEAK_BYTES = 32;
constexpr int TWEAK_NIBBLES = TWEAK_BYTES * 2;
constexpr int KEY_BYTES = 160;
constexpr int KEY_NIBBLES = TWEAK_BYTES * 2;
constexpr int RA = 3;   
constexpr int RB = 5;  


constexpr uint8_t HW2[256] = {0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
                             1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
                             1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
                             0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
                             1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
                             0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
                             0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
                             1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
                             1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
                             0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
                             0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
                             1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
                             0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
                             1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
                             1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
                             0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0};

constexpr uint8_t SBOX[16] = {
    0x1,0x0,0x9,0x3,
    0x8,0x5,0xe,0x7,
    0x4,0x2,0xC,0xB,
    0xA,0xF,0x6,0xD
};


constexpr uint8_t INV_SBOX[16] = {
    0x1,0x0,0x9,0x3,
    0x8,0x5,0xe,0x7,
    0x4,0x2,0xC,0xB,
    0xA,0xF,0x6,0xD
};

constexpr uint8_t M_MATRIX[4][4] = {
    {0, 1, 1, 1},
    {1, 0, 1, 1},
    {1, 1, 0, 1},
    {1, 1, 1, 0}
};

constexpr uint8_t PBOX[STATE_NIBBLES] = {
    5, 12, 4, 1, 17, 9, 10, 16, 28, 14, 21, 22, 11, 27, 8, 13,
    2, 25, 18, 3, 30, 6, 19, 20, 0, 23, 24, 31, 7, 15, 29, 26
};

constexpr uint8_t ROUND_CONST[(RA+RB)][STATE_BYTES] = {
0x44, 0x73, 0x70, 0x03, 0x2e, 0x8a, 0x19, 0x13, 0xd3, 0x08, 0xa3, 0x85, 0x88, 0x6a, 0x3f, 0x24, 
0x89, 0x6c, 0x4e, 0xec, 0x98, 0xfa, 0x2e, 0x08, 0xd0, 0x31, 0x9f, 0x29, 0x22, 0x38, 0x09, 0xa4, 
0x6c, 0x0c, 0xe9, 0x34, 0xcf, 0x66, 0x54, 0xbe, 0x77, 0x13, 0xd0, 0x38, 0xe6, 0x21, 0x28, 0x45, 
0x17, 0x09, 0x47, 0xb5, 0xb5, 0xd5, 0x84, 0x3f, 0xdd, 0x50, 0x7c, 0xc9, 0xb7, 0x29, 0xac, 0xc0, 
0xac, 0xb5, 0xdf, 0x98, 0xa6, 0x0b, 0x31, 0xd1, 0x1b, 0xfb, 0x79, 0x89, 0xd9, 0xd5, 0x16, 0x92, 
0x96, 0x7e, 0x26, 0x6a, 0xed, 0xaf, 0xe1, 0xb8, 0xb7, 0xdf, 0x1a, 0xd0, 0xdb, 0x72, 0xfd, 0x2f,
0xf7, 0x6c, 0x91, 0xb3, 0x47, 0x99, 0xa1, 0x24, 0x99, 0x7f, 0x2c, 0xf1, 0x45, 0x90, 0x7c, 0xba, 
0x69, 0x4e, 0x57, 0x71, 0xd8, 0x20, 0x69, 0x63, 0x16, 0xfc, 0x8e, 0x85, 0xe2, 0xf2, 0x01, 0x08 
};
constexpr uint8_t ROUND_CONST_Prime[(RA+RB)][STATE_BYTES] = {
0x58, 0xb6, 0x8e, 0x72, 0x8f, 0x74, 0x95, 0x0d, 0x7e, 0x3d, 0x93, 0xf4, 0xa3, 0xfe, 0x58, 0xa4, 
0xb5, 0x59, 0x5a, 0xc2, 0x1d, 0xa4, 0x54, 0x7b, 0xee, 0x4a, 0x15, 0x82, 0x58, 0xcd, 0x8b, 0x71, 
0xf0, 0x85, 0x60, 0x28, 0x23, 0xb0, 0xd1, 0xc5, 0x13, 0x60, 0xf2, 0x2a, 0x39, 0xd5, 0x30, 0x9c, 
0x0e, 0x18, 0x3a, 0x60, 0xb0, 0xdc, 0x79, 0x8e, 0xef, 0x38, 0xdb, 0xb8, 0x18, 0x79, 0x41, 0xca, 
0x27, 0x4b, 0x31, 0xbd, 0xc1, 0x77, 0x15, 0xd7, 0x3e, 0x8a, 0x1e, 0xb0, 0x8b, 0x0e, 0x9e, 0x6c, 
0x94, 0xab, 0x55, 0xaa, 0xf3, 0x25, 0x55, 0xe6, 0x60, 0x5c, 0x60, 0x55, 0xda, 0x2f, 0xaf, 0x78,
0xb6, 0x10, 0xab, 0x2a, 0x6a, 0x39, 0xca, 0x55, 0x40, 0x14, 0xe8, 0x63, 0x62, 0x98, 0x48, 0x57, 
0x93, 0xe9, 0x72, 0x7c, 0xaf, 0x86, 0x54, 0xa1, 0xce, 0xe8, 0x41, 0x11, 0x34, 0x5c, 0xcc, 0xb4
}; 



void SubBytes(uint8_t state[STATE_BYTES]) {
    for (int i = 0; i < STATE_BYTES; i++) {
        uint8_t hi = SBOX[state[i] >> 4];
        uint8_t lo = SBOX[state[i] & 0xF];
        state[i] = (hi << 4) | lo;
    }
}

void MixColumns(uint8_t state[STATE_BYTES]) {
    for (int col = 0; col < STATE_NIBBLES / 4; col++) {
        uint8_t coldata[4];
        for (int r = 0; r < 4; r++) {
            int idx = col + r * (STATE_NIBBLES / 4);
            int byteIndex = idx / 2;
            bool highNibble = (idx % 2 == 1);
            uint8_t nibble = highNibble ? (state[byteIndex] >> 4) : (state[byteIndex] & 0xF);
            coldata[r] = nibble;
        }
        uint8_t result[4] = {0};
        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 4; c++) {
                result[r] ^= M_MATRIX[r][c]? coldata[c] : 0;
            }
        }
        for (int r = 0; r < 4; r++) {
            int idx = col + r * (STATE_NIBBLES / 4);
            int byteIndex = idx / 2;
            bool highNibble = (idx % 2 == 1);
            if (highNibble) state[byteIndex] = (result[r] << 4) | (state[byteIndex] & 0xF);
            else state[byteIndex] = (state[byteIndex] & 0xF0) | result[r];
        }
    }
}

void AddRoundKey(uint8_t state[STATE_BYTES], const uint8_t roundKey[STATE_BYTES]) {
    for (int i = 0; i < STATE_BYTES; i++) {
        state[i] ^= roundKey[i];
    }
}

void AddRoundConstant(uint8_t state[STATE_BYTES], const uint8_t constant[STATE_BYTES]) {
    for (int i = 0; i < STATE_BYTES; i++) {
        state[i] ^= constant[i];
    }
}

void Permutation(uint8_t state[STATE_BYTES]) {
    uint8_t temp[STATE_NIBBLES];
    for (int i = 0; i < STATE_NIBBLES; i++) {
        int byteIndex = i / 2;
        bool highNibble = (i % 2 == 1);
        temp[i] = highNibble ? (state[byteIndex] >> 4) : (state[byteIndex] & 0xF);
    }
    uint8_t permuted[STATE_NIBBLES];
    for (int i = 0; i < STATE_NIBBLES; i++) {
        permuted[i] = temp[PBOX[i]];
    }
    for (int i = 0; i < STATE_BYTES; i++) {
        state[i] = (permuted[2*i+1] << 4) | permuted[2*i];
    }
}

void InvPermutation(uint8_t state[STATE_BYTES]) {
    uint8_t temp[STATE_NIBBLES];
    for (int i = 0; i < STATE_NIBBLES; i++) {
        int byteIndex = i / 2;
        bool highNibble = (i % 2 == 1);
        temp[i] = highNibble ? (state[byteIndex] >> 4) : (state[byteIndex] & 0xF);
    }
    uint8_t permuted[STATE_NIBBLES];
    for (int i = 0; i < STATE_NIBBLES; i++) {
        permuted[PBOX[i]] = temp[i];
    }
    for (int i = 0; i < STATE_BYTES; i++) {
        state[i] = (permuted[2*i+1] << 4) | permuted[2*i];
    }
}

void Hash(uint8_t key[STATE_BYTES+TWEAK_BYTES], uint8_t t[TWEAK_BYTES], uint8_t h[STATE_BYTES]) //the key listed in the most (STATE_BYTES+TWEAK_BYTES)*8-1 bits
{    
    for (int i = STATE_BYTES - 1; i >= 0; i--) {
        h[STATE_BYTES - 1 - i] = 0;
        for (int l = 0; l < 8; l++) {
            uint8_t temp[TWEAK_BYTES] = {0};
            for (int j = 0; j < TWEAK_BYTES; j++) {
                temp[TWEAK_BYTES - 1 - j] = (key[TWEAK_BYTES + i - j] << l) ^ (key[TWEAK_BYTES + i - j - 1] >> (8 - l));
            }
            uint8_t p = 0;
            for (int j = 0; j < TWEAK_BYTES; j++) {
                p ^= (t[j] & temp[j]);
            }
            h[STATE_BYTES - 1 - i] ^= (HW2[p] << (l));
        }
    }

}

void GenerateRoundKey(uint8_t masterKey[KEY_BYTES], uint8_t t[TWEAK_BYTES], uint8_t rk[RA+RB][STATE_BYTES], uint8_t w[2][STATE_BYTES], uint8_t h[2][STATE_BYTES]) {
    uint8_t keyPrime[KEY_BYTES] = {0};
    for (int i = 0; i < KEY_BYTES; i++) {
        for (int j = 0; j < 8; j++) {
            keyPrime[i] ^= ((masterKey[(11 * (8 * i + j) % (KEY_BYTES * 8)) / 8] >> ((11 * (8 * i + j) % (KEY_BYTES * 8)) % 8)) & 1) << j;
        }
    }
    for (int i = 0; i < STATE_BYTES; i++) {
        w[0][i] = masterKey[i];
        w[1][i] = masterKey[i+STATE_BYTES];
        for (int j = 0; j < RA+RB; j++) {
            rk[j][i] = masterKey[i + (j+2)*STATE_BYTES];
        }   
    }
    uint8_t hk[2][STATE_BYTES + TWEAK_BYTES];
    for (int i = STATE_BYTES + TWEAK_BYTES - 1; i >= 0; i--){
        if (i > 0) {
            hk[0][i] = (keyPrime[i] << 1) ^ (keyPrime[i - 1] >> 7);;
            hk[1][i] = (keyPrime[i + STATE_BYTES + TWEAK_BYTES] << 2) ^ (keyPrime[i + STATE_BYTES + TWEAK_BYTES - 1] >> 6);
        } else {
            hk[0][i] = keyPrime[i] << 1;
            hk[1][i] = ((keyPrime[i + STATE_BYTES + TWEAK_BYTES] << 2) ^ (keyPrime[i + STATE_BYTES + TWEAK_BYTES - 1] >> 6)) & 0xfe;
        }
    }
    /*cout <<  "hk:" << endl;
    for (int i = 0; i < STATE_BYTES + TWEAK_BYTES; i++) {
        cout << hex << setw(2) << setfill('0') << (int) hk[0][i] << " ";
    }
    cout << endl;
    for (int i = 0; i < STATE_BYTES + TWEAK_BYTES; i++) {
        cout << hex << setw(2) << setfill('0') << (int) hk[1][i] << " ";
    }
    cout << endl;*/
    Hash(hk[0], t, h[0]);
    Hash(hk[1], t, h[1]);
}

void Whitening(uint8_t state[STATE_BYTES], uint8_t w[STATE_BYTES]) {
    for (int i = 0; i < STATE_BYTES; i++)
        state[i] ^= w[i];
}


void Encrypt(uint8_t state[STATE_BYTES], uint8_t rk[RA+RB][STATE_BYTES],
    uint8_t w[2][STATE_BYTES], uint8_t h[2][STATE_BYTES]) {
    Whitening(state, w[0]);
    for (int i = 0; i < STATE_BYTES; i++) cout << hex << setw(2) << setfill('0') << (int)state[i] << " ";
    cout << endl;
    for (int r = 0; r < RA; r++) {
        SubBytes(state);
        MixColumns(state);
        AddRoundKey(state, rk[r]);
        AddRoundConstant(state, ROUND_CONST[r]);
        Permutation(state);
    }
    SubBytes(state);
    MixColumns(state);
    AddRoundKey(state, h[0]);
    Permutation(state);
    for (int r = 0; r < RB; r++) {
        SubBytes(state);
        MixColumns(state);
        AddRoundKey(state, rk[r+RA]);
        AddRoundConstant(state, ROUND_CONST[r+RA]);
        Permutation(state);
    }

    uint8_t h_xor[STATE_BYTES];
    for (int i = 0; i < STATE_BYTES; i++) h_xor[i] = h[0][i] ^ h[1][i];
    SubBytes(state);
    MixColumns(state);
    AddRoundKey(state, h_xor);
    SubBytes(state);

    for (int r = 0; r < RB; r++) {
        InvPermutation(state);
        AddRoundConstant(state, ROUND_CONST_Prime[r]);
        AddRoundKey(state, rk[r]);
        MixColumns(state);
        SubBytes(state);
    }
    InvPermutation(state);
    AddRoundKey(state, h[1]);
    MixColumns(state);
    SubBytes(state);
    for (int r = 0; r < RA; r++) {
        InvPermutation(state);
        AddRoundConstant(state, ROUND_CONST_Prime[r + RB]);
        AddRoundKey(state, rk[r+RB]);
        MixColumns(state);
        SubBytes(state);
    }
    Whitening(state, w[1]);
}

void Decrypt(uint8_t state[STATE_BYTES], uint8_t rk[RA+RB][STATE_BYTES],
    uint8_t w[2][STATE_BYTES], uint8_t h[2][STATE_BYTES]) {
    Whitening(state, w[1]);
    //for (int i = 0; i < STATE_BYTES; i++) cout << hex << setw(2) << setfill('0') << (int)state[i] << " ";
    //cout << endl;
    for (int r = 0; r < RA; r++) {
        SubBytes(state);
        MixColumns(state);
        AddRoundKey(state, rk[RA+RB-r-1]);
        AddRoundConstant(state, ROUND_CONST_Prime[RA+RB-r-1]);
        Permutation(state);
        //for (int i = 0; i < STATE_BYTES; i++) cout << hex << setw(2) << setfill('0') << (int)state[i] << " ";
        //cout << endl;
    }
    SubBytes(state);
    MixColumns(state);
    AddRoundKey(state, h[1]);
    Permutation(state);
    for (int r = 0; r < RB; r++) {
        SubBytes(state);
        MixColumns(state);
        AddRoundKey(state, rk[RB-r-1]);
        AddRoundConstant(state, ROUND_CONST_Prime[RB-r-1]);
        Permutation(state);
    }

    uint8_t h_xor[STATE_BYTES];
    for (int i = 0; i < STATE_BYTES; i++) h_xor[i] = h[0][i] ^ h[1][i];
    SubBytes(state);
    AddRoundKey(state, h_xor);
    MixColumns(state);
    SubBytes(state);

    for (int r = 0; r < RB; r++) {
        InvPermutation(state);
        AddRoundConstant(state, ROUND_CONST[RA+RB-r-1]);
        AddRoundKey(state, rk[RA+RB-r-1]);
        MixColumns(state);
        SubBytes(state);
    }
    InvPermutation(state);
    AddRoundKey(state, h[0]);
    MixColumns(state);
    SubBytes(state);
    for (int r = 0; r < RA; r++) {
        InvPermutation(state);
        AddRoundConstant(state, ROUND_CONST[RA - r - 1]);
        AddRoundKey(state, rk[RA-r-1]);
        MixColumns(state);
        SubBytes(state);
    }
    Whitening(state, w[0]);
}

int main() {
    uint8_t state[STATE_BYTES] = {0};
    uint8_t masterKey[KEY_BYTES]={0x8d, 0xd2, 0x48, 0x4a, 0x5a, 0x3f, 0xd7, 0x13, 0x08, 0xc1, 0x2f, 0x3b, 0xe4, 0x3b, 0x96, 0xd6, 0xc2, 0x95, 0x63, 0x6a, 0x02, 0x17, 0x2c, 0x75, 0xda, 0x3e, 0x89, 0x96, 0x4c, 0x2a, 0x96, 0x28, 0xe4, 0x66, 0xa6, 0x88, 0xd8, 0x02, 0xa1, 0xd6, 0xbb, 0x74, 0x43, 0x90, 0x04, 0xf1, 0xb8, 0x21, 0xe4, 0x64, 0xa6, 0x34, 0xe2, 0x02, 0xa1, 0xd6, 0xbb, 0x74, 0x43, 0x90, 0x04, 0xf1, 0xb8, 0x21, 0xb6, 0x39, 0xae, 0x51, 0x02, 0xbd, 0xa0, 0x57, 0xca, 0x42, 0xa3, 0xe8, 0x0f, 0x30, 0x40, 0xcc, 0x45, 0x0a, 0xc3, 0x69, 0xb7, 0x0c, 0xd7, 0x74, 0xe1, 0xba, 0x78, 0xd5, 0x48, 0xc6, 0x61, 0x29, 0x0f, 0xf8, 0x17, 0x05, 0x35, 0x66, 0x2b, 0x5e, 0xa1, 0x7f, 0x8e, 0xb3, 0x21, 0x90, 0x77, 0x97, 0x06, 0x8a, 0xc7, 0x78, 0xad, 0xe0, 0x22, 0x30, 0xce, 0x52, 0x0b, 0x4a, 0x87, 0x87, 0xd3, 0x6d, 0x21, 0x0d, 0x6f, 0xde, 0xf9, 0x2d, 0x2c, 0x76, 0x7e, 0xc5, 0x02, 0xf3, 0xc6, 0x1d, 0x7c, 0xe0, 0x43, 0x69, 0x24, 0x3a, 0xc3, 0xde, 0xd7, 0xd1, 0xe4, 0x67, 0xa4, 0x88, 0xd8, 0x02, 0xa1, 0xd6};
    uint8_t t[TWEAK_BYTES]={0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01, 0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01, 0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01, 0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01};
    uint8_t rk[RA+RB][STATE_BYTES];
    uint8_t w[2][STATE_BYTES];
    uint8_t h[2][STATE_BYTES];

    GenerateRoundKey(masterKey, t, rk, w, h);
    cout << "Plaintext: ";
    for (auto b : state) cout << hex << setw(2) << setfill('0') << (int)b << " ";
    cout << "\n";



    Encrypt(state, rk, w, h);

    cout << "Ciphertext: ";
    for (auto b : state) cout << hex << setw(2) << setfill('0') << (int)b << " ";
    cout << "\n";

    /*for (int i = 0; i < 2; i++) {
        for (auto b : h[i]) cout << hex << setw(2) << setfill('0') << (int)b << " ";
        cout << endl;
    }*/

    Decrypt(state, rk, w, h);

    cout << "Recovered Plaintext: ";
    for (auto b : state) cout << hex << setw(2) << setfill('0') << (int)b << " ";
    cout << "\n";
}
