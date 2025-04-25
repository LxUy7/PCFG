#include "md5.h"
#include <iomanip>
#include <cassert> // Use <cassert> for assert
#include <chrono>
#include <vector>
#include <stdexcept> // For throwing exceptions

using namespace std;
using namespace std::chrono; // Use std::chrono

// --- Original StringProcess function (unchanged) ---
/**
 * StringProcess: 将单个输入字符串转换成MD5计算所需的消息数组
 * @param input 输入
 * @param[out] n_byte 用于给调用者传递额外的返回值，即最终Byte数组的长度
 * @return Byte消息数组 (调用者需要负责释放内存)
 */
Byte *StringProcess(const string& input, size_t *n_byte) // Use size_t for length, const& for input
{
    // 将输入的字符串转换为Byte为单位的数组
    // 注意：直接使用 c_str() 的指针并不安全，如果 input 被销毁，指针会失效。
    // 最好是复制数据。但这里原始代码就是这样用的，暂时保留。
    // 更安全的做法：vector<Byte> data(input.begin(), input.end());
    const Byte *blocks = (const Byte *)input.c_str(); // Input data pointer
    size_t length = input.length();                  // Use size_t

    // 计算原始消息长度（以比特为单位）
    uint64_t bitLength = (uint64_t)length * 8; // Use uint64_t for bit length

    // paddingBits: 原始消息需要的padding长度（以bit为单位）
    // 计算填充使得总长度 (bits) 模 512 等于 448
    size_t paddingBits = 0;
    size_t remainder = bitLength % 512;

    if (remainder >= 448) {
        paddingBits = (512 - remainder) + 448;
    } else {
        paddingBits = 448 - remainder;
    }
    // 至少填充 1 bit (0x80 byte)
    if (paddingBits == 0) {
        paddingBits = 512;
    }


    // 原始消息需要的padding长度（以Byte为单位）
    size_t paddingBytes = paddingBits / 8;
    assert(paddingBytes > 0); // Must have at least one padding byte (0x80)

    // 创建最终的字节数组
    // length + paddingBytes + 8:
    // 1. length为原始消息的长度（Bytes）
    // 2. paddingBytes为原始消息需要的padding长度（Bytes） (包括0x80)
    // 3. 额外附加64bits的原始消息长度，即8个bytes
    size_t paddedLength = length + paddingBytes + 8;
    Byte *paddedMessage = new Byte[paddedLength];

    // 复制原始消息
    memcpy(paddedMessage, blocks, length);

    // 添加填充字节。填充时，第一位为1，后面的所有位均为0。
    // 所以第一个byte是0x80
    paddedMessage[length] = 0x80;
    // 填充0字节 (注意 paddingBytes 包含了 0x80 这一字节)
    if (paddingBytes > 1) {
         memset(paddedMessage + length + 1, 0, paddingBytes - 1);
    }

    // 添加消息长度（64比特，小端格式）
    for (int i = 0; i < 8; ++i)
    {
        // 这里的length是原始消息的长度 (bytes) * 8 = bitLength
        paddedMessage[length + paddingBytes + i] = (Byte)((bitLength >> (i * 8)) & 0xFF);
    }

    // 验证长度是否满足要求。此时长度应当是512bit (64 bytes) 的倍数
    assert((paddedLength * 8) % 512 == 0);
    assert(paddedLength % 64 == 0);

    // 在填充+添加长度之后，消息被分为n_blocks个512bit的部分
    *n_byte = paddedLength;
    return paddedMessage;
}


// --- Original MD5Hash function (unchanged, for scalar/fallback) ---
/**
 * MD5Hash: 将单个输入字符串转换成MD5 (标量版本)
 * @param input 输入
 * @param[out] state 用于给调用者传递额外的返回值，即最终的缓冲区 (长度为4)，也就是MD5的结果
 */
void MD5Hash(const string& input, bit32 *state) // Use const& for input
{
    Byte *paddedMessage = nullptr; // Initialize pointer
    size_t messageLength = 0;      // Use size_t, initialize

    // --- 准备消息 ---
    // The loop `for (int i = 0; i < 1; i += 1)` in the original code
    // was redundant. Calling StringProcess once is enough.
    paddedMessage = StringProcess(input, &messageLength);
    if (!paddedMessage) {
        // Handle allocation failure or error from StringProcess if it could return null
        throw std::runtime_error("Failed to process input string for MD5.");
    }

    assert((messageLength % 64) == 0); // 确保长度是64字节的倍数
    size_t n_blocks = messageLength / 64;

    // --- 初始化状态 ---
    state[0] = 0x67452301;
    state[1] = 0xefcdab89;
    state[2] = 0x98badcfe;
    state[3] = 0x10325476;

    // --- 逐block地更新state ---
    for (size_t i = 0; i < n_blocks; ++i) // Use size_t for loop counter
    {
        bit32 x[16]; // 当前块的16个32位字 (小端)
        const Byte* blockPtr = paddedMessage + i * 64;

        // --- 从块中解码16个字 (小端) ---
        for (int j = 0; j < 16; ++j)
        {
            const Byte* p = blockPtr + j * 4;
            x[j] = (bit32)p[0] | ((bit32)p[1] << 8) | ((bit32)p[2] << 16) | ((bit32)p[3] << 24);
        }

        // --- 保存当前状态 ---
        bit32 a = state[0];
        bit32 b = state[1];
        bit32 c = state[2];
        bit32 d = state[3];

        // --- 执行4轮共64步计算 (使用原始宏) ---
        // **注意：原始代码在这里使用了 #define 宏，这些宏现在已被注释掉或未包含**
        // **要使此原始函数工作，你需要取消注释 md5.h 中的原始 #define F,G,H,I,ROTATELEFT,FF,GG,HH,II**
        // **或者提供这些宏的定义。以下代码假定这些宏是可用的。**

        /* Round 1 */
        #define F(x, y, z) (((x) & (y)) | ((~x) & (z)))
        #define ROTATELEFT(num, n) (((num) << (n)) | ((num) >> (32-(n))))
        #define FF(a, b, c, d, x, s, ac) { \
        (a) += F ((b), (c), (d)) + (x) + ac; \
        (a) = ROTATELEFT ((a), (s)); \
        (a) += (b); \
        }
        FF(a, b, c, d, x[0], s11, 0xd76aa478);
        FF(d, a, b, c, x[1], s12, 0xe8c7b756);
        FF(c, d, a, b, x[2], s13, 0x242070db);
        FF(b, c, d, a, x[3], s14, 0xc1bdceee);
        FF(a, b, c, d, x[4], s11, 0xf57c0faf);
        FF(d, a, b, c, x[5], s12, 0x4787c62a);
        FF(c, d, a, b, x[6], s13, 0xa8304613);
        FF(b, c, d, a, x[7], s14, 0xfd469501);
        FF(a, b, c, d, x[8], s11, 0x698098d8);
        FF(d, a, b, c, x[9], s12, 0x8b44f7af);
        FF(c, d, a, b, x[10], s13, 0xffff5bb1);
        FF(b, c, d, a, x[11], s14, 0x895cd7be);
        FF(a, b, c, d, x[12], s11, 0x6b901122);
        FF(d, a, b, c, x[13], s12, 0xfd987193);
        FF(c, d, a, b, x[14], s13, 0xa679438e);
        FF(b, c, d, a, x[15], s14, 0x49b40821);

        /* Round 2 */
         #define G(x, y, z) (((x) & (z)) | ((y) & (~z)))
         #define GG(a, b, c, d, x, s, ac) { \
        (a) += G ((b), (c), (d)) + (x) + ac; \
        (a) = ROTATELEFT ((a), (s)); \
        (a) += (b); \
        }
        GG(a, b, c, d, x[1], s21, 0xf61e2562);
        GG(d, a, b, c, x[6], s22, 0xc040b340);
        GG(c, d, a, b, x[11], s23, 0x265e5a51);
        GG(b, c, d, a, x[0], s24, 0xe9b6c7aa);
        GG(a, b, c, d, x[5], s21, 0xd62f105d);
        GG(d, a, b, c, x[10], s22, 0x2441453);
        GG(c, d, a, b, x[15], s23, 0xd8a1e681);
        GG(b, c, d, a, x[4], s24, 0xe7d3fbc8);
        GG(a, b, c, d, x[9], s21, 0x21e1cde6);
        GG(d, a, b, c, x[14], s22, 0xc33707d6);
        GG(c, d, a, b, x[3], s23, 0xf4d50d87);
        GG(b, c, d, a, x[8], s24, 0x455a14ed);
        GG(a, b, c, d, x[13], s21, 0xa9e3e905);
        GG(d, a, b, c, x[2], s22, 0xfcefa3f8);
        GG(c, d, a, b, x[7], s23, 0x676f02d9);
        GG(b, c, d, a, x[12], s24, 0x8d2a4c8a);

        /* Round 3 */
        #define H(x, y, z) ((x) ^ (y) ^ (z))
        #define HH(a, b, c, d, x, s, ac) { \
        (a) += H ((b), (c), (d)) + (x) + ac; \
        (a) = ROTATELEFT ((a), (s)); \
        (a) += (b); \
        }
        HH(a, b, c, d, x[5], s31, 0xfffa3942);
        HH(d, a, b, c, x[8], s32, 0x8771f681);
        HH(c, d, a, b, x[11], s33, 0x6d9d6122);
        HH(b, c, d, a, x[14], s34, 0xfde5380c);
        HH(a, b, c, d, x[1], s31, 0xa4beea44);
        HH(d, a, b, c, x[4], s32, 0x4bdecfa9);
        HH(c, d, a, b, x[7], s33, 0xf6bb4b60);
        HH(b, c, d, a, x[10], s34, 0xbebfbc70);
        HH(a, b, c, d, x[13], s31, 0x289b7ec6);
        HH(d, a, b, c, x[0], s32, 0xeaa127fa);
        HH(c, d, a, b, x[3], s33, 0xd4ef3085);
        HH(b, c, d, a, x[6], s34, 0x4881d05);
        HH(a, b, c, d, x[9], s31, 0xd9d4d039);
        HH(d, a, b, c, x[12], s32, 0xe6db99e5);
        HH(c, d, a, b, x[15], s33, 0x1fa27cf8);
        HH(b, c, d, a, x[2], s34, 0xc4ac5665);

        /* Round 4 */
         #define I(x, y, z) ((y) ^ ((x) | (~z)))
         #define II(a, b, c, d, x, s, ac) { \
        (a) += I ((b), (c), (d)) + (x) + ac; \
        (a) = ROTATELEFT ((a), (s)); \
        (a) += (b); \
        }
        II(a, b, c, d, x[0], s41, 0xf4292244);
        II(d, a, b, c, x[7], s42, 0x432aff97);
        II(c, d, a, b, x[14], s43, 0xab9423a7);
        II(b, c, d, a, x[5], s44, 0xfc93a039);
        II(a, b, c, d, x[12], s41, 0x655b59c3);
        II(d, a, b, c, x[3], s42, 0x8f0ccc92);
        II(c, d, a, b, x[10], s43, 0xffeff47d);
        II(b, c, d, a, x[1], s44, 0x85845dd1);
        II(a, b, c, d, x[8], s41, 0x6fa87e4f);
        II(d, a, b, c, x[15], s42, 0xfe2ce6e0);
        II(c, d, a, b, x[6], s43, 0xa3014314);
        II(b, c, d, a, x[13], s44, 0x4e0811a1);
        II(a, b, c, d, x[4], s41, 0xf7537e82);
        II(d, a, b, c, x[11], s42, 0xbd3af235);
        II(c, d, a, b, x[2], s43, 0x2ad7d2bb);
        II(b, c, d, a, x[9], s44, 0xeb86d391);
        // --- 添加结果到状态 ---
        state[0] += a;
        state[1] += b;
        state[2] += c;
        state[3] += d;
    }

    // --- Endian Swap (字节序转换) ---
    // 这个操作在原始代码中有，但通常MD5的内部状态是小端，
    // 字节序转换通常只在最终输出为十六进制字符串时进行。
    // 这里保留原始行为。
    for (int i = 0; i < 4; i++)
    {
        uint32_t value = state[i];
        state[i] = ((value & 0xff) << 24) |       // 将最低字节移到最高位
                   ((value & 0xff00) << 8) |      // 将次低字节左移
                   ((value & 0xff0000) >> 8) |    // 将次高字节右移
                   ((value & 0xff000000) >> 24);  // 将最高字节移到最低位
    }

    cout<<"使用了原本的"<<endl;
    // 释放动态分配的内存
    delete[] paddedMessage;
    // delete[] messageLength; // messageLength was not allocated with new[]
}


// --- NEON SIMD Version for processing 4 inputs ---

// 辅助函数：从块中加载并解码4个输入的第j个字到NEON向量
// 使用临时数组 + vld1q_u32 的方式进行 gather
static inline bit32x4 load_x_vec(const Byte* p0, const Byte* p1, const Byte* p2, const Byte* p3, int j) {
    // 确保数组对齐以获得最佳性能
    alignas(16) uint32_t vals[4];

    const Byte* ptr0 = p0 + j * 4;
    const Byte* ptr1 = p1 + j * 4;
    const Byte* ptr2 = p2 + j * 4;
    const Byte* ptr3 = p3 + j * 4;

    // 手动小端解码
    vals[0] = (uint32_t)ptr0[0] | ((uint32_t)ptr0[1] << 8) | ((uint32_t)ptr0[2] << 16) | ((uint32_t)ptr0[3] << 24);
    vals[1] = (uint32_t)ptr1[0] | ((uint32_t)ptr1[1] << 8) | ((uint32_t)ptr1[2] << 16) | ((uint32_t)ptr1[3] << 24);
    vals[2] = (uint32_t)ptr2[0] | ((uint32_t)ptr2[1] << 8) | ((uint32_t)ptr2[2] << 16) | ((uint32_t)ptr2[3] << 24);
    vals[3] = (uint32_t)ptr3[0] | ((uint32_t)ptr3[1] << 8) | ((uint32_t)ptr3[2] << 16) | ((uint32_t)ptr3[3] << 24);

    return vld1q_u32(vals); // 从对齐的数组加载到向量
}

/**
 * MD5Hash_NEON_4x: 使用NEON指令并行计算4个输入字符串的MD5
 * @param inputs 包含4个输入字符串的vector (必须正好是4个)
 * @param[out] outputs 指向一个4x4的bit32数组，用于存储4个MD5结果 (每个结果4个bit32)
 * @note **重要假设**: 所有4个输入经过填充后具有相同的块数。
 */
void MD5Hash_NEON_4x(const std::vector<string>& inputs, bit32 (*outputs)[4])
{
    if (inputs.size() != 4) {
        throw std::invalid_argument("MD5Hash_NEON_4x requires exactly 4 input strings.");
    }

    Byte* paddedMessages[4] = {nullptr, nullptr, nullptr, nullptr};
    size_t messageLengths[4] = {0};
    size_t n_blocks = 0;

    // --- 准备4个消息并检查长度一致性 ---
    try {
        for (int i = 0; i < 4; ++i) {
            paddedMessages[i] = StringProcess(inputs[i], &messageLengths[i]);
            if (!paddedMessages[i]) {
                 throw std::runtime_error("Failed to process input string for MD5.");
            }
            if (i == 0) {
                assert((messageLengths[0] % 64) == 0);
                n_blocks = messageLengths[0] / 64;
            } else {
                // 检查填充后的长度是否一致
                if (messageLengths[i] != messageLengths[0]) {
                    throw std::runtime_error("MD5Hash_NEON_4x requires all inputs to have the same padded length.");
                }
            }
        }
    } catch (...) {
        // 清理已分配的内存
        for(int i=0; i<4; ++i) delete[] paddedMessages[i];
        throw; // 重新抛出异常
    }


    // --- 初始化4组状态向量 ---
    bit32x4 v_a = vdupq_n_u32(0x67452301);
    bit32x4 v_b = vdupq_n_u32(0xefcdab89);
    bit32x4 v_c = vdupq_n_u32(0x98badcfe);
    bit32x4 v_d = vdupq_n_u32(0x10325476);

    // --- 预先创建常量向量 ---
    // Round 1 Constants
    const bit32x4 K0  = vdupq_n_u32(0xd76aa478); const bit32x4 K1  = vdupq_n_u32(0xe8c7b756);
    const bit32x4 K2  = vdupq_n_u32(0x242070db); const bit32x4 K3  = vdupq_n_u32(0xc1bdceee);
    const bit32x4 K4  = vdupq_n_u32(0xf57c0faf); const bit32x4 K5  = vdupq_n_u32(0x4787c62a);
    const bit32x4 K6  = vdupq_n_u32(0xa8304613); const bit32x4 K7  = vdupq_n_u32(0xfd469501);
    const bit32x4 K8  = vdupq_n_u32(0x698098d8); const bit32x4 K9  = vdupq_n_u32(0x8b44f7af);
    const bit32x4 K10 = vdupq_n_u32(0xffff5bb1); const bit32x4 K11 = vdupq_n_u32(0x895cd7be);
    const bit32x4 K12 = vdupq_n_u32(0x6b901122); const bit32x4 K13 = vdupq_n_u32(0xfd987193);
    const bit32x4 K14 = vdupq_n_u32(0xa679438e); const bit32x4 K15 = vdupq_n_u32(0x49b40821);
    // Round 2 Constants
    const bit32x4 K16 = vdupq_n_u32(0xf61e2562); const bit32x4 K17 = vdupq_n_u32(0xc040b340);
    const bit32x4 K18 = vdupq_n_u32(0x265e5a51); const bit32x4 K19 = vdupq_n_u32(0xe9b6c7aa);
    const bit32x4 K20 = vdupq_n_u32(0xd62f105d); const bit32x4 K21 = vdupq_n_u32(0x02441453); // Note: leading zero for literal
    const bit32x4 K22 = vdupq_n_u32(0xd8a1e681); const bit32x4 K23 = vdupq_n_u32(0xe7d3fbc8);
    const bit32x4 K24 = vdupq_n_u32(0x21e1cde6); const bit32x4 K25 = vdupq_n_u32(0xc33707d6);
    const bit32x4 K26 = vdupq_n_u32(0xf4d50d87); const bit32x4 K27 = vdupq_n_u32(0x455a14ed);
    const bit32x4 K28 = vdupq_n_u32(0xa9e3e905); const bit32x4 K29 = vdupq_n_u32(0xfcefa3f8);
    const bit32x4 K30 = vdupq_n_u32(0x676f02d9); const bit32x4 K31 = vdupq_n_u32(0x8d2a4c8a);
    // Round 3 Constants
    const bit32x4 K32 = vdupq_n_u32(0xfffa3942); const bit32x4 K33 = vdupq_n_u32(0x8771f681);
    const bit32x4 K34 = vdupq_n_u32(0x6d9d6122); const bit32x4 K35 = vdupq_n_u32(0xfde5380c);
    const bit32x4 K36 = vdupq_n_u32(0xa4beea44); const bit32x4 K37 = vdupq_n_u32(0x4bdecfa9);
    const bit32x4 K38 = vdupq_n_u32(0xf6bb4b60); const bit32x4 K39 = vdupq_n_u32(0xbebfbc70);
    const bit32x4 K40 = vdupq_n_u32(0x289b7ec6); const bit32x4 K41 = vdupq_n_u32(0xeaa127fa);
    const bit32x4 K42 = vdupq_n_u32(0xd4ef3085); const bit32x4 K43 = vdupq_n_u32(0x04881d05); // Note: leading zero
    const bit32x4 K44 = vdupq_n_u32(0xd9d4d039); const bit32x4 K45 = vdupq_n_u32(0xe6db99e5);
    const bit32x4 K46 = vdupq_n_u32(0x1fa27cf8); const bit32x4 K47 = vdupq_n_u32(0xc4ac5665);
    // Round 4 Constants
    const bit32x4 K48 = vdupq_n_u32(0xf4292244); const bit32x4 K49 = vdupq_n_u32(0x432aff97);
    const bit32x4 K50 = vdupq_n_u32(0xab9423a7); const bit32x4 K51 = vdupq_n_u32(0xfc93a039);
    const bit32x4 K52 = vdupq_n_u32(0x655b59c3); const bit32x4 K53 = vdupq_n_u32(0x8f0ccc92);
    const bit32x4 K54 = vdupq_n_u32(0xffeff47d); const bit32x4 K55 = vdupq_n_u32(0x85845dd1);
    const bit32x4 K56 = vdupq_n_u32(0x6fa87e4f); const bit32x4 K57 = vdupq_n_u32(0xfe2ce6e0);
    const bit32x4 K58 = vdupq_n_u32(0xa3014314); const bit32x4 K59 = vdupq_n_u32(0x4e0811a1);
    const bit32x4 K60 = vdupq_n_u32(0xf7537e82); const bit32x4 K61 = vdupq_n_u32(0xbd3af235);
    const bit32x4 K62 = vdupq_n_u32(0x2ad7d2bb); const bit32x4 K63 = vdupq_n_u32(0xeb86d391);


    // --- 逐block地更新状态向量 ---
    for (size_t i = 0; i < n_blocks; ++i)
    {
        // --- 加载当前块的数据到16个向量 ---
        bit32x4 x_vec[16];
        const Byte* p0 = paddedMessages[0] + i * 64;
        const Byte* p1 = paddedMessages[1] + i * 64;
        const Byte* p2 = paddedMessages[2] + i * 64;
        const Byte* p3 = paddedMessages[3] + i * 64;

        for (int j = 0; j < 16; ++j) {
            x_vec[j] = load_x_vec(p0, p1, p2, p3, j);
        }

        // --- 保存当前状态向量 ---
        bit32x4 v_aa = v_a;
        bit32x4 v_bb = v_b;
        bit32x4 v_cc = v_c;
        bit32x4 v_dd = v_d;

        // --- 执行4轮共64步计算 (使用NEON内联函数) ---
        /* Round 1 */
        FF_neon<s11>(v_a, v_b, v_c, v_d, x_vec[ 0], K0);
        FF_neon<s12>(v_d, v_a, v_b, v_c, x_vec[ 1], K1);
        FF_neon<s13>(v_c, v_d, v_a, v_b, x_vec[ 2], K2);
        FF_neon<s14>(v_b, v_c, v_d, v_a, x_vec[ 3], K3);
        FF_neon<s11>(v_a, v_b, v_c, v_d, x_vec[ 4], K4);
        FF_neon<s12>(v_d, v_a, v_b, v_c, x_vec[ 5], K5);
        FF_neon<s13>(v_c, v_d, v_a, v_b, x_vec[ 6], K6);
        FF_neon<s14>(v_b, v_c, v_d, v_a, x_vec[ 7], K7);
        FF_neon<s11>(v_a, v_b, v_c, v_d, x_vec[ 8], K8);
        FF_neon<s12>(v_d, v_a, v_b, v_c, x_vec[ 9], K9);
        FF_neon<s13>(v_c, v_d, v_a, v_b, x_vec[10], K10);
        FF_neon<s14>(v_b, v_c, v_d, v_a, x_vec[11], K11);
        FF_neon<s11>(v_a, v_b, v_c, v_d, x_vec[12], K12);
        FF_neon<s12>(v_d, v_a, v_b, v_c, x_vec[13], K13);
        FF_neon<s13>(v_c, v_d, v_a, v_b, x_vec[14], K14);
        FF_neon<s14>(v_b, v_c, v_d, v_a, x_vec[15], K15);

        /* Round 2 */
        GG_neon<s21>(v_a, v_b, v_c, v_d, x_vec[ 1], K16);
        GG_neon<s22>(v_d, v_a, v_b, v_c, x_vec[ 6], K17);
        GG_neon<s23>(v_c, v_d, v_a, v_b, x_vec[11], K18);
        GG_neon<s24>(v_b, v_c, v_d, v_a, x_vec[ 0], K19);
        GG_neon<s21>(v_a, v_b, v_c, v_d, x_vec[ 5], K20);
        GG_neon<s22>(v_d, v_a, v_b, v_c, x_vec[10], K21);
        GG_neon<s23>(v_c, v_d, v_a, v_b, x_vec[15], K22);
        GG_neon<s24>(v_b, v_c, v_d, v_a, x_vec[ 4], K23);
        GG_neon<s21>(v_a, v_b, v_c, v_d, x_vec[ 9], K24);
        GG_neon<s22>(v_d, v_a, v_b, v_c, x_vec[14], K25);
        GG_neon<s23>(v_c, v_d, v_a, v_b, x_vec[ 3], K26);
        GG_neon<s24>(v_b, v_c, v_d, v_a, x_vec[ 8], K27);
        GG_neon<s21>(v_a, v_b, v_c, v_d, x_vec[13], K28);
        GG_neon<s22>(v_d, v_a, v_b, v_c, x_vec[ 2], K29);
        GG_neon<s23>(v_c, v_d, v_a, v_b, x_vec[ 7], K30);
        GG_neon<s24>(v_b, v_c, v_d, v_a, x_vec[12], K31);

        /* Round 3 */
        HH_neon<s31>(v_a, v_b, v_c, v_d, x_vec[ 5], K32);
        HH_neon<s32>(v_d, v_a, v_b, v_c, x_vec[ 8], K33);
        HH_neon<s33>(v_c, v_d, v_a, v_b, x_vec[11], K34);
        HH_neon<s34>(v_b, v_c, v_d, v_a, x_vec[14], K35);
        HH_neon<s31>(v_a, v_b, v_c, v_d, x_vec[ 1], K36);
        HH_neon<s32>(v_d, v_a, v_b, v_c, x_vec[ 4], K37);
        HH_neon<s33>(v_c, v_d, v_a, v_b, x_vec[ 7], K38);
        HH_neon<s34>(v_b, v_c, v_d, v_a, x_vec[10], K39);
        HH_neon<s31>(v_a, v_b, v_c, v_d, x_vec[13], K40);
        HH_neon<s32>(v_d, v_a, v_b, v_c, x_vec[ 0], K41);
        HH_neon<s33>(v_c, v_d, v_a, v_b, x_vec[ 3], K42);
        HH_neon<s34>(v_b, v_c, v_d, v_a, x_vec[ 6], K43);
        HH_neon<s31>(v_a, v_b, v_c, v_d, x_vec[ 9], K44);
        HH_neon<s32>(v_d, v_a, v_b, v_c, x_vec[12], K45);
        HH_neon<s33>(v_c, v_d, v_a, v_b, x_vec[15], K46);
        HH_neon<s34>(v_b, v_c, v_d, v_a, x_vec[ 2], K47);

        /* Round 4 */
        II_neon<s41>(v_a, v_b, v_c, v_d, x_vec[ 0], K48);
        II_neon<s42>(v_d, v_a, v_b, v_c, x_vec[ 7], K49);
        II_neon<s43>(v_c, v_d, v_a, v_b, x_vec[14], K50);
        II_neon<s44>(v_b, v_c, v_d, v_a, x_vec[ 5], K51);
        II_neon<s41>(v_a, v_b, v_c, v_d, x_vec[12], K52);
        II_neon<s42>(v_d, v_a, v_b, v_c, x_vec[ 3], K53);
        II_neon<s43>(v_c, v_d, v_a, v_b, x_vec[10], K54);
        II_neon<s44>(v_b, v_c, v_d, v_a, x_vec[ 1], K55);
        II_neon<s41>(v_a, v_b, v_c, v_d, x_vec[ 8], K56);
        II_neon<s42>(v_d, v_a, v_b, v_c, x_vec[15], K57);
        II_neon<s43>(v_c, v_d, v_a, v_b, x_vec[ 6], K58);
        II_neon<s44>(v_b, v_c, v_d, v_a, x_vec[13], K59);
        II_neon<s41>(v_a, v_b, v_c, v_d, x_vec[ 4], K60);
        II_neon<s42>(v_d, v_a, v_b, v_c, x_vec[11], K61);
        II_neon<s43>(v_c, v_d, v_a, v_b, x_vec[ 2], K62);
        II_neon<s44>(v_b, v_c, v_d, v_a, x_vec[ 9], K63);


        // --- 添加结果到状态向量 ---
        v_a = vaddq_u32(v_a, v_aa);
        v_b = vaddq_u32(v_b, v_bb);
        v_c = vaddq_u32(v_c, v_cc);
        v_d = vaddq_u32(v_d, v_dd);
        // cout<<"使用了neon"<<endl;
    }

    // --- 提取最终结果到输出数组 ---
    outputs[0][0] = vgetq_lane_u32(v_a, 0); outputs[0][1] = vgetq_lane_u32(v_b, 0); outputs[0][2] = vgetq_lane_u32(v_c, 0); outputs[0][3] = vgetq_lane_u32(v_d, 0);
    outputs[1][0] = vgetq_lane_u32(v_a, 1); outputs[1][1] = vgetq_lane_u32(v_b, 1); outputs[1][2] = vgetq_lane_u32(v_c, 1); outputs[1][3] = vgetq_lane_u32(v_d, 1);
    outputs[2][0] = vgetq_lane_u32(v_a, 2); outputs[2][1] = vgetq_lane_u32(v_b, 2); outputs[2][2] = vgetq_lane_u32(v_c, 2); outputs[2][3] = vgetq_lane_u32(v_d, 2);
    outputs[3][0] = vgetq_lane_u32(v_a, 3); outputs[3][1] = vgetq_lane_u32(v_b, 3); outputs[3][2] = vgetq_lane_u32(v_c, 3); outputs[3][3] = vgetq_lane_u32(v_d, 3);

    // --- Endian Swap (字节序转换) ---
    // 应用于每个输出结果，以匹配原始 scalar 函数的行为
    for(int i=0; i<4; ++i) {
        for (int j = 0; j < 4; j++)
        {
            uint32_t value = outputs[i][j];
            outputs[i][j] = ((value & 0xff) << 24) |
                            ((value & 0xff00) << 8) |
                            ((value & 0xff0000) >> 8) |
                            ((value & 0xff000000) >> 24);
        }
    }

    // --- 清理内存 ---
    for(int i=0; i<4; ++i) {
        delete[] paddedMessages[i];
    }
}