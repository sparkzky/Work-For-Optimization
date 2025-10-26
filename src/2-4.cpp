#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

/**
 * @brief 根据输入的种子生成一个伪随机浮点数 (程序示例2.39)
 */
float rand_float(float s) { return 4.0f * s * (1.0f - s); }

/**
 * @brief 生成两个 N*N 的浮点矩阵 a 和 b (程序示例2.39)
 */
void matrix_gen(float* a, float* b, int N, float seed) {
    float s = seed;
    long long size = (long long)N * N;
    for (long long i = 0; i < size; ++i) {
        s = rand_float(s);
        a[i] = s;
        s = rand_float(s);
        b[i] = s;
    }
}

float calculate_trace(const float* c, int N) {
    float trace = 0.0f;
    for (int i = 0; i < N; ++i) {
        trace += c[i * N + i];
    }
    return trace;
}

template <typename MultiplyFunc>
std::chrono::duration<double>
run_matrix_multiply_test(const std::string& name, int N, float seed, MultiplyFunc multiply_func) {
    std::cout << name << " N=" << N << " seed=" << seed << std::endl;

    std::vector<float> a((long long)N * N);
    std::vector<float> b((long long)N * N);
    std::vector<float> c((long long)N * N);

    matrix_gen(a.data(), b.data(), N, seed);

    auto start = std::chrono::high_resolution_clock::now();
    multiply_func(a.data(), b.data(), c.data(), N);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;
    float trace = calculate_trace(c.data(), N);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Trace: " << trace << std::endl;
    std::cout << "计算时间(s): " << duration.count() << std::endl;

    return duration;
}

// =========================================================================
//                      作业任务 (2): 基准矩阵乘法
// =========================================================================

/**
 * @brief 基准矩阵乘法实现 (ijk 顺序) (程序示例2.40)
 */
void matrix_multiply_baseline(float* a, float* b, float* c, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            float sum = 0.0f;
            for (int k = 0; k < N; ++k) {
                sum += a[i * N + k] * b[k * N + j];
            }
            c[i * N + j] = sum;
        }
    }
}

std::chrono::duration<double> basic_multiply(int n, float seed) {
    if (n <= 0 || seed <= 0.0f || seed >= 1.0f) {
        std::cerr << "错误: 参数无效。请确保 N > 0 且 0 < seed < 1。" << std::endl;
        return std::chrono::duration<double>::zero();
    }

    return run_matrix_multiply_test("Matrix_mul", n, seed, matrix_multiply_baseline);
}

// 基本矩阵乘法测试
void test_basic_multiply(float seed = 0.12345f) {
    const int sizes[] = {512, 1024, 2048, 4096};
    for (int n : sizes) {
        basic_multiply(n, seed);
    }
}

// =========================================================================
//                   作业任务 (5): 分块矩阵乘法
// =========================================================================

/**
 * @brief 使用分块技术进行矩阵乘法
 */
void matrix_multiply_blocked(float* a, float* b, float* c, int N, int m) {
    // 预先清零 C 矩阵
    for (long long i = 0; i < (long long)N * N; ++i) {
        c[i] = 0.0f;
    }

    for (int i0 = 0; i0 < N; i0 += m) {
        for (int j0 = 0; j0 < N; j0 += m) {
            for (int k0 = 0; k0 < N; k0 += m) {
                int i_limit = std::min(i0 + m, N);
                for (int i = i0; i < i_limit; ++i) {
                    int j_limit = std::min(j0 + m, N);
                    for (int j = j0; j < j_limit; ++j) {
                        int k_limit = std::min(k0 + m, N);
                        float sum_val = c[i * N + j];
                        for (int k = k0; k < k_limit; ++k) {
                            sum_val += a[i * N + k] * b[k * N + j];
                        }
                        c[i * N + j] = sum_val;
                    }
                }
            }
        }
    }
}

std::chrono::duration<double> blocked_multiply(int block_size, int N = 4096, float seed = 0.12345f) {
    return run_matrix_multiply_test(
        "Blocked_multiply (block=" + std::to_string(block_size) + ")",
        N,
        seed,
        [block_size](float* a, float* b, float* c, int N) { matrix_multiply_blocked(a, b, c, N, block_size); });
}

// =========================================================================
//                                主函数
// =========================================================================

void print_usage(const char* prog_name) {
    std::cerr << "用法: " << prog_name << " <task_id> [options]" << std::endl;
    std::cerr << "任务列表:" << std::endl;
    std::cerr << "  1: 回答理论问题" << std::endl;
    std::cerr << "  2: 运行基准矩阵乘法 (N=512, 1024, 2048, 4096)" << std::endl;
    std::cerr << "  5: 设计矩阵分块程序 (N=4096, m=m0/2, m0, 2m0)" << std::endl;
    std::cerr << "     用法: " << prog_name << " 5 <m0>" << std::endl;
    std::cerr << "  6: 在最优分块尺寸下运行，并计算加速比和MFLOPS" << std::endl;
    std::cerr << "     用法: " << prog_name << " 6 <optimal_m>" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    int task_id = std::atoi(argv[1]);
    const float seed = 0.12345f;

    switch (task_id) {
    case 1: {
        std::cout << "========== 作业问题 (1), (3), (4) 理论解答 ==========" << std::endl;
        std::cout << "\n--- 问题(1): 两个N阶矩阵乘法的计算复杂度是多少？所需浮点计算有多少次？---\n";
        std::cout << "时间复杂度: O(N^3)\n";
        std::cout << "浮点计算次数: 对于结果矩阵C的每个元素C[i][j]，需要进行N次乘法和N-1次加法。";
        std::cout << "为了简化，通常将一次乘加视为一次核心计算，即N次乘加。";
        std::cout
            << "总共有 N*N 个元素，所以总计算量约为 N*N*N 次乘加，即 2*N^3 次浮点运算（N^3次乘法 + N^3次加法）。\n";
        std::cout << "例如，N=512时，浮点运算次数为 2 * 512^3 = 268,435,456 次。\n";

        std::cout << "\n--- 问题(3): 运行时间变化是否符合理论复杂度的预期？ ---\n";
        std::cout << "是的。矩阵乘法的时间复杂度是 O(N^3)。当N加倍时（例如从512到1024），";
        std::cout << "理论上运行时间应该增加到原来的 (2N)^3 / N^3 = 8 倍。";
        std::cout
            << "实际测试中，由于缓存效应和内存带宽的限制，实际增加的倍数可能在6到8倍之间，但总体趋势是符合O(N^3)的。\n";

        std::cout << "\n--- 问题(4): 硬件平台的一级数据Cache容量多少？分块时m取多少合适？ ---\n";
        std::cout << "请根据你自己的CPU信息填写，例如：\n";
        std::cout << "我的CPU是 Intel Core i5-8500，其L1数据缓存为每个核心32KB。\n";
        std::cout << "采用分块方法时，我们希望三个m*m的子块（一个来自A，一个来自B，一个用于C）能完全放入L1缓存。";
        std::cout << "所需内存为 3 * m * m * sizeof(float) = 12 * m^2 字节。\n";
        std::cout << "为使其小于32KB (32768字节)，12 * m^2 < 32768，解得 m < 52.1。\n";
        std::cout << "因此，选择一个接近52且便于计算（如2的幂次）的值是合适的，例如 m=32 或 m=48。\n";
        break;
    }

    case 2: {
        std::cout << "========== 作业任务 (2): 测试基准矩阵乘法运行时间 ==========" << std::endl;
        test_basic_multiply();
        break;
    }

    case 5: {
        if (argc != 3) {
            std::cerr << "错误: 任务5需要提供一个m0值。\n";
            print_usage(argv[0]);
            return 1;
        }
        int m0 = std::atoi(argv[2]);
        const int N = 4096;
        std::vector<int> ms = {m0 / 2, m0, m0 * 2};

        std::cout << "========== 作业任务 (5): 测试不同分块尺寸 (m0=" << m0 << ") ==========" << std::endl;

        double min_time = -1.0;
        int best_m = 0;

        for (int m : ms) {
            if (m == 0)
                continue;
            auto duration = blocked_multiply(m);

            if (min_time < 0 || duration.count() < min_time) {
                min_time = duration.count();
                best_m = m;
            }
        }
        std::cout << "\n-----------------------------------------------------" << std::endl;
        std::cout << "结论: 时间最短的分块方法是 m = " << best_m << ", 时间为 " << min_time << "s." << std::endl;
        break;
    }

    case 6: {
        if (argc != 3) {
            std::cerr << "错误: 任务6需要提供一个最优m值。\n";
            print_usage(argv[0]);
            return 1;
        }
        int optimal_m = std::atoi(argv[2]);
        const int N = 4096;

        std::cout << "========== 作业任务 (6): 在最优分块尺寸下性能分析 (m=" << optimal_m
                  << ") ==========" << std::endl;

        // 运行基准版本
        std::cout << "\n--- 运行基准版本 (N=4096) ---" << std::endl;
        auto time_base = basic_multiply(4096, seed);

        // 运行最优分块版本
        std::cout << "\n--- 运行最优分块版本 (N=4096, m=" << optimal_m << ") ---" << std::endl;
        auto time_opt = blocked_multiply(optimal_m);

        // 计算并输出结果
        double speedup = time_base.count() / time_opt.count();
        // 浮点运算次数 = 2 * N^3
        double flop = 2.0 * pow(N, 3);
        // MFLOPS = (浮点运算次数 / 1,000,000) / 时间
        double mflops = (flop / 1e6) / time_opt.count();

        std::cout << "\n-----------------------------------------------------" << std::endl;
        std::cout << "性能分析结果:" << std::endl;
        std::cout << "加速比 (相比基准程序): " << std::fixed << std::setprecision(2) << speedup << " 倍" << std::endl;
        std::cout << "每秒钟完成的单精度浮点计算次数: " << mflops << " MFLOPS" << std::endl;

        break;
    }

    default:
        std::cerr << "错误: 无效的任务ID。" << std::endl;
        print_usage(argv[0]);
        return 1;
    }

    return 0;
}