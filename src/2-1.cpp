#include "../include/2-1.h"

#include <chrono>  // 用于高精度计时
#include <iomanip> // 用于设置输出精度
#include <iostream>

// ======================================================================
// 测例 1: 标准 Horner 算法 (无优化)
// ======================================================================

double calculate_s_standard() {
    double sum = 0.0;
    for (long long k = 0; k < N; ++k) {
        sum += f_horner(k * h);
    }
    return sum * h;
}

// ======================================================================
// 测例 2: 4次循环展开
// ======================================================================

double calculate_s_unrolled_x4() {
    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    long long limit = N - (N % 4);

    for (long long k = 0; k < limit; k += 4) {
        sum0 += f_horner((k + 0) * h);
        sum1 += f_horner((k + 1) * h);
        sum2 += f_horner((k + 2) * h);
        sum3 += f_horner((k + 3) * h);
    }

    // 清理收尾
    for (long long k = limit; k < N; ++k) {
        sum0 += f_horner(k * h);
    }

    return (sum0 + sum1 + sum2 + sum3) * h;
}

// ======================================================================
// 测例 3: 使用 AVX 指令 (不展开)
// ======================================================================

double calculate_s_avx() {
    // __m256d 是一个可以容纳 4 个 double 的向量类型
    __m256d sum_vec = _mm256_setzero_pd(); // 向量累加器，初始化为 [0.0, 0.0, 0.0, 0.0]

    // 将标量常量广播到向量的所有通道
    const __m256d h_vec = _mm256_set1_pd(h);
    const __m256d C0_v = _mm256_set1_pd(C0);
    const __m256d C1_v = _mm256_set1_pd(C1);
    const __m256d C2_v = _mm256_set1_pd(C2);
    const __m256d C3_v = _mm256_set1_pd(C3);
    // 用于步进 k 的向量 [0, 1, 2, 3]
    const __m256d k_step = _mm256_set_pd(3.0, 2.0, 1.0, 0.0);

    long long limit = N - (N % 4);
    for (long long k = 0; k < limit; k += 4) {
        __m256d k_vec = _mm256_add_pd(_mm256_set1_pd((double)k), k_step);
        __m256d x_vec = _mm256_mul_pd(k_vec, h_vec);
        __m256d f_vec = f_horner_avx(x_vec, C0_v, C1_v, C2_v, C3_v);
        sum_vec = _mm256_add_pd(sum_vec, f_vec);
    }

    // 水平求和: 将向量 [s3, s2, s1, s0] 内的 4 个元素相加
    double sum_partials[4];
    _mm256_storeu_pd(sum_partials, sum_vec);
    double total_sum = sum_partials[0] + sum_partials[1] + sum_partials[2] + sum_partials[3];

    // 清理收尾 (标量计算)
    for (long long k = limit; k < N; ++k) {
        total_sum += f_horner(k * h);
    }

    return total_sum * h;
}

// ======================================================================
// 测例 4: 使用 AVX 指令 + 4次循环展开
// ======================================================================

double calculate_s_avx_unrolled_x4() {
    // 4个向量累加器，打破数据依赖
    __m256d sum_vec0 = _mm256_setzero_pd();
    __m256d sum_vec1 = _mm256_setzero_pd();
    __m256d sum_vec2 = _mm256_setzero_pd();
    __m256d sum_vec3 = _mm256_setzero_pd();

    const __m256d h_vec = _mm256_set1_pd(h);
    const __m256d C0_v = _mm256_set1_pd(C0);
    const __m256d C1_v = _mm256_set1_pd(C1);
    const __m256d C2_v = _mm256_set1_pd(C2);
    const __m256d C3_v = _mm256_set1_pd(C3);
    const __m256d k_step = _mm256_set_pd(3.0, 2.0, 1.0, 0.0);
    const __m256d k_unroll_step = _mm256_set1_pd(4.0); // 每次展开处理4个数

    // 每次迭代处理 4*4 = 16 个元素
    long long limit = N - (N % 16);
    for (long long k = 0; k < limit; k += 16) {
        // 第一次迭代: k to k+3
        __m256d k_base0 = _mm256_set1_pd((double)k);
        __m256d x_vec0 = _mm256_mul_pd(_mm256_add_pd(k_base0, k_step), h_vec);
        sum_vec0 = _mm256_add_pd(sum_vec0, f_horner_avx(x_vec0, C0_v, C1_v, C2_v, C3_v));

        // 第二次迭代: k+4 to k+7
        __m256d k_base1 = _mm256_add_pd(k_base0, k_unroll_step);
        __m256d x_vec1 = _mm256_mul_pd(_mm256_add_pd(k_base1, k_step), h_vec);
        sum_vec1 = _mm256_add_pd(sum_vec1, f_horner_avx(x_vec1, C0_v, C1_v, C2_v, C3_v));

        // 第三次迭代: k+8 to k+11
        __m256d k_base2 = _mm256_add_pd(k_base1, k_unroll_step);
        __m256d x_vec2 = _mm256_mul_pd(_mm256_add_pd(k_base2, k_step), h_vec);
        sum_vec2 = _mm256_add_pd(sum_vec2, f_horner_avx(x_vec2, C0_v, C1_v, C2_v, C3_v));

        // 第四次迭代: k+12 to k+15
        __m256d k_base3 = _mm256_add_pd(k_base2, k_unroll_step);
        __m256d x_vec3 = _mm256_mul_pd(_mm256_add_pd(k_base3, k_step), h_vec);
        sum_vec3 = _mm256_add_pd(sum_vec3, f_horner_avx(x_vec3, C0_v, C1_v, C2_v, C3_v));
    }

    // 合并4个向量累加器
    __m256d total_vec_sum = _mm256_add_pd(_mm256_add_pd(sum_vec0, sum_vec1), _mm256_add_pd(sum_vec2, sum_vec3));

    double sum_partials[4];
    _mm256_storeu_pd(sum_partials, total_vec_sum);
    double total_sum = sum_partials[0] + sum_partials[1] + sum_partials[2] + sum_partials[3];

    // 清理收尾
    for (long long k = limit; k < N; ++k) {
        total_sum += f_horner(k * h);
    }

    return total_sum * h;
}

int main() {
    std::cout << std::fixed << std::setprecision(15);
    std::cout << "实验环境: N = " << N << ", h = " << h << std::endl;

    auto run_test = [](const char* name, double (*func)()) {
        std::cout << "\n-------------------------------------------\n";
        std::cout << "正在运行: " << name << "...\n";
        auto start = std::chrono::high_resolution_clock::now();
        double s = func();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "计算结果 s = " << s << std::endl;
        std::cout << "执行时间: " << duration.count() << " 秒" << std::endl;
        return duration.count();
    };

    double t1 = run_test("测例 1: 标准 Horner 算法", calculate_s_standard);
    double t2 = run_test("测例 2: Horner + 4次循环展开", calculate_s_unrolled_x4);
    double t3 = run_test("测例 3: AVX (不展开)", calculate_s_avx);
    double t4 = run_test("测例 4: AVX + 4次循环展开", calculate_s_avx_unrolled_x4);

    std::cout << "\n-------------------------------------------\n";
    std::cout << "性能对比总结:\n";
    std::cout << " - 循环展开 vs 标准: " << t1 / t2 << " 倍速\n";
    std::cout << " - AVX vs 标准:      " << t1 / t3 << " 倍速\n";
    std::cout << " - AVX+展开 vs 标准: " << t1 / t4 << " 倍速\n";
    std::cout << " - AVX+展开 vs AVX:  " << t3 / t4 << " 倍速\n";

    return 0;
}
