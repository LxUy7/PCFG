#include "PCFG.h"
#include <chrono>
#include <fstream>
#include "md5.h"
#include <iomanip>
using namespace std;
using namespace chrono;

// 编译指令如下：
// g++ correctness.cpp train.cpp guessing.cpp md5.cpp -o test.exe
// arm-linux-gnueabihf-g++ correctness.cpp train.cpp guessing.cpp md5.cpp -o test.exe
// cd /mnt/d/大二下/并行/作业三/PCFG_framework
//arm-linux-gnueabihf-g++ -mfpu=neon -mfloat-abi=hard -march=armv7-a correctness.cpp train.cpp guessing.cpp md5.cpp -o test.exe
// ./test.exe
//qemu-arm -L /usr/arm-linux-gnueabihf/ ./test.exe

// 通过这个函数，你可以验证你实现的SIMD哈希函数的正确性
int main()
{
    bit32 state[4];
    MD5Hash_NEON_4x("bvaisdbjasdkafkasdfnavkjnakdjfejfanjsdnfkajdfkajdfjkwanfdjaknsvjkanbjbjadfajwefajksdfakdnsvjadfasjdvabvaisdbjasdkafkasdfnavkjnakdjfejfanjsdnfkajdfkajdfjkwanfdjaknsvjkanbjbjadfajwefajksdfakdnsvjadfasjdvabvaisdbjasdkafkasdfnavkjnakdjfejfanjsdnfkajdfkajdfjkwanfdjaknsvjkanbjbjadfajwefajksdfakdnsvjadfasjdvabvaisdbjasdkafkasdfnavkjnakdjfejfanjsdnfkajdfkajdfjkwanfdjaknsvjkanbjbjadfajwefajksdfakdnsvjadfasjdva", state);
    for (int i1 = 0; i1 < 4; i1 += 1)
    {
        cout << std::setw(8) << std::setfill('0') << hex << state[i1];
    }
    cout << endl;
}
/*bba46eb8b53cf65d50ca54b2f8afd9db
*/