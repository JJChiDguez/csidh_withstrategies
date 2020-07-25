#ifndef _ADDC_H_
#define _ADDC_H_

// The public Edwards curve with coefficients a and b. E is isomorphic to E / F_p : y^2 = x^3 + x
// a and (a - d) in Montgomery domain represetantion 
static proj E = {
{ 0x767762E5FD1E1599, 0x33C5743A49A0B6F6, 0x68FC0C0364C77443, 0xB9AA1E24F83F56DB, 0x3914101F20520EFB, 0x7B1ED6D95B1542B4, 0x114A8BE928C8828A, 0x3793732BBB24F40},	// a
{ 0xECEEC5CBFA3C2B32, 0x678AE87493416DEC, 0xD1F81806C98EE886, 0x73543C49F07EADB6, 0x7228203E40A41DF7, 0xF63DADB2B62A8568, 0x229517D251910514, 0x6F26E6577649E80}	// a - d
};

// Shortest differential addition chains for each l_i
static uint64_t ADDITION_CHAIN[] = {
   0x0,   0x0,   0x1,   0x1,   0x0,   0x5,   0x4,   0xD,
   0xA,   0x8,   0x3,  0x18,  0x11,   0x1,  0x36,  0x32,
  0x2C,   0x9,  0x24,  0x14,   0x4,  0x6C,   0x0,  0x49,
  0x68,   0xC,  0x60,  0x4A,  0xE1,  0xE8,  0x40,  0x2B,
  0xC5,  0xD4,  0x98,  0xC2,  0xD0,  0xA1,  0xC0,  0x8A,
  0x81,  0x88,  0x82, 0x198,   0x1, 0x191, 0x185, 0x194,
 0x134,   0x0, 0x181, 0x1A2, 0x118, 0x109, 0x12A,  0x51,
 0x141, 0x122, 0x148, 0x144, 0x101, 0x108,  0x14,  0x50,
 0x140, 0x2D8,  0x10, 0x324, 0x231, 0x352,  0xD1, 0x312,
 0x268, 0x612
};

// Length of the shortest differential addition chain for each l_i
static uint8_t ADDITION_CHAIN_LENGTH[] = {
  0,  1,  2,  3,  3,  4,  4,  5,
  5,  5,  6,  6,  6,  6,  7,  7,
  7,  7,  7,  7,  7,  8,  7,  8,
  8,  8,  8,  8,  9,  9,  8,  9,
  9,  9,  9,  9,  9,  9,  9,  9,
  9,  9,  9, 10,  9, 10, 10, 10,
 10,  9, 10, 10, 10, 10, 10, 10,
 10, 10, 10, 10, 10, 10, 10, 10,
 10, 11, 10, 11, 11, 11, 11, 11,
 11, 12
};

// L
static uint32_t L[] = { 
   3,   5,   7,  11,  13,  17,  19,  23,
  29,  31,  37,  41,  43,  47,  53,  59,
  61,  67,  71,  73,  79,  83,  89,  97,
 101, 103, 107, 109, 113, 127, 131, 137,
 139, 149, 151, 157, 163, 167, 173, 179,
 181, 191, 193, 197, 199, 211, 223, 227,
 229, 233, 239, 241, 251, 257, 263, 269,
 271, 277, 281, 283, 293, 307, 311, 313,
 317, 331, 337, 347, 349, 353, 359, 367,
 373, 587
};

static uint16_t BITS_OF_L[] = { 
  2,  3,  3,  4,  4,  5,  5,  5,
  5,  5,  6,  6,  6,  6,  6,  6,
  6,  7,  7,  7,  7,  7,  7,  7,
  7,  7,  7,  7,  7,  7,  8,  8,
  8,  8,  8,  8,  8,  8,  8,  8,
  8,  8,  8,  8,  8,  8,  8,  8,
  8,  8,  8,  8,  8,  9,  9,  9,
  9,  9,  9,  9,  9,  9,  9,  9,
  9,  9,  9,  9,  9,  9,  9,  9,
  9, 10
};

#define BITS_OF_4SQRT_OF_P 258
#define LARGE_L 587
// The l_is are only required for isogeny constructions
#endif /* Addition chains */
