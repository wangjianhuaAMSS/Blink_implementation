#include <iostream>
#include <iomanip>
#include <stdio.h>
using namespace std;

constexpr int STATE_BYTES = 8;
constexpr int STATE_NIBBLES = STATE_BYTES * 2;
constexpr int TWEAK_BYTES = 8;
constexpr int TWEAK_NIBBLES = TWEAK_BYTES * 2;
constexpr int KEY_BYTES = 56;
constexpr int KEY_NIBBLES = TWEAK_BYTES * 2;
constexpr int RA = 2;   
constexpr int RB = 3;  


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
    0, 5, 11, 10, 
    1, 6, 4, 13, 
    2, 12, 9, 15, 
    3, 7, 14, 8
};

constexpr uint8_t ROUND_CONST[(RA+RB)][STATE_BYTES] = {
0x44, 0x73, 0x70, 0x03, 0x2e, 0x8a, 0x19, 0x13, 
0x89, 0x6c, 0x4e, 0xec, 0x98, 0xfa, 0x2e, 0x08, 
0x6c, 0x0c, 0xe9, 0x34, 0xcf, 0x66, 0x54, 0xbe, 
0x17, 0x09, 0x47, 0xb5, 0xb5, 0xd5, 0x84, 0x3f, 
0xac, 0xb5, 0xdf, 0x98, 0xa6, 0x0b, 0x31, 0xd1};
constexpr uint8_t ROUND_CONST_Prime[(RA+RB)][STATE_BYTES] = {
0x58, 0xb6, 0x8e, 0x72, 0x8f, 0x74, 0x95, 0x0d,
0xb5, 0x59, 0x5a, 0xc2, 0x1d, 0xa4, 0x54, 0x7b, 
0xf0, 0x85, 0x60, 0x28, 0x23, 0xb0, 0xd1, 0xc5, 
0x0e, 0x18, 0x3a, 0x60, 0xb0, 0xdc, 0x79, 0x8e, 
0x27, 0x4b, 0x31, 0xbd, 0xc1, 0x77, 0x15, 0xd7
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
    cout <<  "hk:" << endl;
    for (int i = 0; i < STATE_BYTES + TWEAK_BYTES; i++) {
        cout << hex << setw(2) << setfill('0') << (int) hk[0][i] << " ";
    }
    cout << endl;
    for (int i = 0; i < STATE_BYTES + TWEAK_BYTES; i++) {
        cout << hex << setw(2) << setfill('0') << (int) hk[1][i] << " ";
    }
    cout << endl;
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
    //for (int i = 0; i < STATE_BYTES; i++) cout << hex << setw(2) << setfill('0') << (int)state[i] << " ";
    //cout << endl;
    for (int r = 0; r < RA; r++) {
        SubBytes(state);
        MixColumns(state);
        AddRoundKey(state, rk[r]);
        AddRoundConstant(state, ROUND_CONST[r]);
        Permutation(state);
        //for (int i = 0; i < STATE_BYTES; i++) cout << hex << setw(2) << setfill('0') << (int)state[i] << " ";
        //cout << endl;
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
    uint8_t masterKey[KEY_BYTES]={0xa1, 0x7f, 0x8e, 0xb3, 0x21, 0x90, 0x77, 0x97, 0x06, 0x8a, 0xc7, 0x78, 0xad, 0xe0, 0x22, 0x30, 0xce, 0x52, 0x0b, 0x4a, 0x87, 0x87, 0xd3, 0x6d, 0x21, 0x0d, 0x6f, 0xde, 0xf9, 0x2d, 0x2c, 0x76, 0x7e, 0xc5, 0x02, 0xf3, 0xc6, 0x1d, 0x7c, 0xe0, 0x43, 0x69, 0x24, 0x3a, 0xc3, 0xde, 0xd7, 0xd1, 0xe4, 0x67, 0xa4, 0x88, 0xd8, 0x02, 0xa1, 0xd6};
    uint8_t t[TWEAK_BYTES]={0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01};
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

    cout << "Recoverd Plaintext: ";
    for (auto b : state) cout << hex << setw(2) << setfill('0') << (int)b << " ";
    cout << "\n";
}
