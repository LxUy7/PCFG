#include "PCFG.h"
#include <chrono>
#include <fstream>
#include "md5.h"
#include <iomanip>
using namespace std;
using namespace chrono;

// 编译指令如下
// g++ main.cpp train.cpp guessing.cpp md5.cpp -o test.exe
// g++ main.cpp train.cpp guessing.cpp md5.cpp -o test.exe -O1
// g++ main.cpp train.cpp guessing.cpp md5.cpp -o test.exe -O2
// arm-linux-gnueabihf-g++ -mfpu=neon -mfloat-abi=hard -march=armv7-a main.cpp train.cpp guessing.cpp md5.cpp -o test.exe
// qemu-arm -L /usr/arm-linux-gnueabihf/ ./test.exe
// cd /mnt/d/大二下/并行/作业三/PCFG_framework
// arm-linux-gnueabihf-g++ -mfpu=neon -mfloat-abi=hard -march=armv7-a -fopenmp main.cpp train.cpp guessing.cpp md5.cpp -o test.exe

int main()
{
    double time_hash = 0;  // 用于MD5哈希的时间
    double time_guess = 0; // 哈希和猜测的总时长
    double time_train = 0; // 模型训练的总时长
    PriorityQueue q;
    auto start_train = system_clock::now();
    q.m.train("./input/Rockyou-singleLined-full.txt");
    q.m.order();
    auto end_train = system_clock::now();
    auto duration_train = duration_cast<microseconds>(end_train - start_train);
    time_train = double(duration_train.count()) * microseconds::period::num / microseconds::period::den;

    q.init();
    cout << "here" << endl;
    int curr_num = 0;
    auto start = system_clock::now();
    // 由于需要定期清空内存，我们在这里记录已生成的猜测总数
    int history = 0;
    // std::ofstream a("./output/results.txt");
    while (!q.priority.empty())
    {
        q.PopNext();
        q.total_guesses = q.guesses.size();
        if (q.total_guesses - curr_num >= 100000)
        {
            cout << std::dec <<"Guesses generated: " << history + q.total_guesses << endl;
            curr_num = q.total_guesses;

            // 在此处更改实验生成的猜测上限
            int generate_n = 10000000;
            if (history + q.total_guesses >= 10000000)
            {
                auto end = system_clock::now();
                auto duration = duration_cast<microseconds>(end - start);
                time_guess = double(duration.count()) * microseconds::period::num / microseconds::period::den;
                cout << "Guess time:" << time_guess - time_hash << " seconds" << endl;
                cout << "Hash time:" << time_hash << " seconds" << endl;
                cout << "Train time:" << time_train << " seconds" << endl;
                break;
            }
        }

        // 为了避免内存超限，我们在q.guesses中口令达到一定数目时，将其中的所有口令取出并且进行哈希
        // 然后，q.guesses将会被清空。为了有效记录已经生成的口令总数，维护一个history变量来进行记录
        if (curr_num > 1000000)
        {
            auto start_hash = system_clock::now();
            bit32 state[4][4]; // 4个结果，每个结果包含4个32位的MD5哈希

            // 确保每次都处理4个密码，使用 MD5Hash_NEON_4x 进行并行处理
            vector<string> guesses_batch;

            // 处理不足4个的情况
            for (int i = 0; i < 4 && !q.guesses.empty(); i++) {
                guesses_batch.push_back(q.guesses.back());
                q.guesses.pop_back();
            }

            if (guesses_batch.size() == 4) {
                // 使用并行化版本的MD5进行哈希计算
                MD5Hash_NEON_4x(guesses_batch, state);
                
                // 处理输出（或者将其写入文件）
                for (int i = 0; i < 4; i++) {
                    // 输出哈希值，或者其他操作
                    cout << "Password guess " << i << ": " << guesses_batch[i] << " -> MD5: ";
                    for (int j = 0; j < 4; j++) {
                        cout << std::setw(8) << std::setfill('0') << hex << state[i][j];
                    }
                    cout << endl;
                }
            } else {
                // 对于不足4个的情况，调用标量版本的MD5进行处理
                for (const string& pw : guesses_batch) {
                    MD5Hash(pw, state[0]); // 将哈希结果存入 state[0] 或者其他索引
                    // 输出结果
                    cout << "Password guess: " << pw << " -> MD5: ";
                    for (int j = 0; j < 4; j++) {
                        cout << std::setw(8) << std::setfill('0') << hex << state[0][j];
                    }
                    cout << endl;
                }
            }

            // 在这里对哈希所需的总时长进行计算
            auto end_hash = system_clock::now();
            auto duration = duration_cast<microseconds>(end_hash - start_hash);
            time_hash += double(duration.count()) * microseconds::period::num / microseconds::period::den;

            // 记录已经生成的口令总数
            history += curr_num;
            curr_num = 0;
            q.guesses.clear();
        }
    }

    return 0;
}
