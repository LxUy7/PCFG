#ifndef MD5_H
#define MD5_H

#include <iostream>
#include <string>
#include <cstring>
#include <vector> // 包含 vector 以便在 MD5Hash 声明中使用
#include <arm_neon.h> // 包含 NEON 头文件

using namespace std;

// 定义了Byte，便于使用
typedef unsigned char Byte;
// 定义了32比特
typedef unsigned int bit32;

// ---- NEON 相关类型定义 ----
// 定义 128 位 NEON 向量类型，包含 4 个 32 位无符号整数
typedef uint32x4_t bit32x4;

// ---- MD5 核心参数 (保持不变) ----
// MD5的一系列参数。参数是固定的
#define s11 7
#define s12 12
#define s13 17
#define s14 22
#define s21 5
#define s22 9
#define s23 14
#define s24 20
#define s31 4
#define s32 11
#define s33 16
#define s34 23
#define s41 6
#define s42 10
#define s43 15
#define s44 21

// ---- MD5 基础函数 (宏定义，用于标量计算) ----
// 这些是原始的、非 SIMD 的函数宏定义
// F(x, y, z) = (x & y) | (~x & z)
#define F(x, y, z) (((x) & (y)) | ((~x) & (z)))
// G(x, y, z) = (x & z) | (y & ~z)
#define G(x, y, z) (((x) & (z)) | ((y) & (~z)))
// H(x, y, z) = x ^ y ^ z
#define H(x, y, z) ((x) ^ (y) ^ (z))
// I(x, y, z) = y ^ (x | ~z)
#define I(x, y, z) ((y) ^ ((x) | (~z)))

// 原始左循环移位 (宏定义，用于标量计算)
#define ROTATELEFT(num, n) (((num) << (n)) | ((num) >> (32-(n)))) 

// 原始 MD5 轮函数 (宏定义，用于标量计算)
// 注意：这些宏直接修改第一个参数 a
#define FF(a, b, c, d, x, s, ac) { \
  (a) += F ((b), (c), (d)) + (x) + ac; \
  (a) = ROTATELEFT ((a), (s)); \
  (a) += (b); \
}

#define GG(a, b, c, d, x, s, ac) { \
  (a) += G ((b), (c), (d)) + (x) + ac; \
  (a) = ROTATELEFT ((a), (s)); \
  (a) += (b); \
}
#define HH(a, b, c, d, x, s, ac) { \
  (a) += H ((b), (c), (d)) + (x) + ac; \
  (a) = ROTATELEFT ((a), (s)); \
  (a) += (b); \
}
#define II(a, b, c, d, x, s, ac) { \
  (a) += I ((b), (c), (d)) + (x) + ac; \
  (a) = ROTATELEFT ((a), (s)); \
  (a) += (b); \
}

// ---- NEON SIMD 并行化函数 ----

// NEON 版本的 F 函数: (x & y) | (~x & z)
static inline bit32x4 F_neon(bit32x4 x, bit32x4 y, bit32x4 z) {
    return vbslq_u32(x, y, z);
}

// NEON 版本的 G 函数: (x & z) | (y & ~z)
static inline bit32x4 G_neon(bit32x4 x, bit32x4 y, bit32x4 z) {
    return vbslq_u32(z, x, y);
}

// NEON 版本的 H 函数: x ^ y ^ z
static inline bit32x4 H_neon(bit32x4 x, bit32x4 y, bit32x4 z) {
    return veorq_u32(veorq_u32(x, y), z);
}

// NEON 版本的 I 函数: y ^ (x | ~z)
static inline bit32x4 I_neon(bit32x4 x, bit32x4 y, bit32x4 z) {
    return veorq_u32(y, vorrq_u32(x, vmvnq_u32(z)));
}

// NEON 版本的左循环移位函数
template<int n>
static inline bit32x4 ROTATELEFT_neon(bit32x4 num) {
    static_assert(n >= 0 && n < 32, "Rotate amount must be between 0 and 31");
    if (n == 0) return num; // 移位 0 位等于不移位
    return vorrq_u32(vshlq_n_u32(num, n), vshrq_n_u32(num, 32 - n));
}

// NEON 版本的 FF 轮函数
template<int s>
static inline void FF_neon(bit32x4 &a, bit32x4 b, bit32x4 c, bit32x4 d, bit32x4 x, bit32x4 ac) {
    a = vaddq_u32(a, F_neon(b, c, d));
    a = vaddq_u32(a, x);
    a = vaddq_u32(a, ac);
    a = ROTATELEFT_neon<s>(a);
    a = vaddq_u32(a, b);
}

// NEON 版本的 GG 轮函数
template<int s>
static inline void GG_neon(bit32x4 &a, bit32x4 b, bit32x4 c, bit32x4 d, bit32x4 x, bit32x4 ac) {
    a = vaddq_u32(a, G_neon(b, c, d));
    a = vaddq_u32(a, x);
    a = vaddq_u32(a, ac);
    a = ROTATELEFT_neon<s>(a);
    a = vaddq_u32(a, b);
}

// NEON 版本的 HH 轮函数
template<int s>
static inline void HH_neon(bit32x4 &a, bit32x4 b, bit32x4 c, bit32x4 d, bit32x4 x, bit32x4 ac) {
    a = vaddq_u32(a, H_neon(b, c, d));
    a = vaddq_u32(a, x);
    a = vaddq_u32(a, ac);
    a = ROTATELEFT_neon<s>(a);
    a = vaddq_u32(a, b);
}

// NEON 版本的 II 轮函数
template<int s>
static inline void II_neon(bit32x4 &a, bit32x4 b, bit32x4 c, bit32x4 d, bit32x4 x, bit32x4 ac) {
    a = vaddq_u32(a, I_neon(b, c, d));
    a = vaddq_u32(a, x);
    a = vaddq_u32(a, ac);
    a = ROTATELEFT_neon<s>(a);
    a = vaddq_u32(a, b);
}

// ---- 原始 MD5 哈希函数声明 ----
// 计算单个输入的 MD5 哈希值
// 输入：input (字符串), state (长度为4的bit32数组，用于存储结果 A, B, C, D)
void MD5Hash(const string& input, bit32 *state);

// ---- NEON 版本的 MD5 哈希函数声明 ----
// 计算 4 个口令的 MD5 哈希值（NEON 版本）
void MD5Hash_NEON_4x(const std::vector<string>& inputs, bit32 (*states)[4]);

#endif // MD5_H
