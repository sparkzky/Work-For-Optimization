#include <cmath>       // 用于 M_PI
#include <immintrin.h> // 引入 Intel SIMD 指令集头文件 (AVX)


// --- 全局常量定义 ---
const long long N = 10000000000;
const double PI = M_PI;
const double h = (PI / 2.0) / N;

// --- 多项式 f(x) 的系数 ---
const double C3 = -1.0 / 5040.0; // -1/7!
const double C2 = 1.0 / 120.0;   //  1/5!
const double C1 = -1.0 / 6.0;    // -1/3!
const double C0 = 1.0;

// ======================================================================
// 测例 1: 标准 Horner 算法 (无优化)
// ======================================================================

/**
 * @brief 使用霍纳法则计算函数 f(x) 的值 (标量版本)
 * @param x 输入的自变量
 * @return f(x) 的计算结果
 */
inline double f_horner(double x) {
    double x2 = x * x;
    // 应用霍纳法则: C0 + x2*(C1 + x2*(C2 + x2*C3))
    // 使用 FMA (Fused Multiply-Add) 可以提高精度和速度，现代编译器在-O3下通常会自动优化
    double poly_val = std::fma(std::fma(std::fma(C3, x2, C2), x2, C1), x2, C0);
    return x * poly_val;
}

double calculate_s_standard() ;

// ======================================================================
// 测例 2: 4次循环展开
// ======================================================================

double calculate_s_unrolled_x4();

// ======================================================================
// 测例 3: 使用 AVX 指令 (不展开)
// ======================================================================

/**
 * @brief 使用 AVX 指令集并行计算 4 个 f(x) 的值
 * @param x_vec 一个包含 4 个 double 的 AVX 向量 [x, x+h, x+2h, x+3h]
 * @param C0_v, C1_v, C2_v, C3_v 系数的 AVX 向量版本
 * @return 包含 4 个 f(x) 结果的 AVX 向量
 */
inline __m256d f_horner_avx(__m256d x_vec, __m256d C0_v, __m256d C1_v, __m256d C2_v, __m256d C3_v) {
    __m256d x2_vec = _mm256_mul_pd(x_vec, x_vec); // x*x

    // 使用 FMA 指令: _mm256_fmadd_pd(a, b, c) -> (a*b)+c
    __m256d term3 = _mm256_fmadd_pd(C3_v, x2_vec, C2_v);  // C3*x2 + C2
    __m256d term2 = _mm256_fmadd_pd(term3, x2_vec, C1_v); // (C3*x2+C2)*x2 + C1
    __m256d term1 = _mm256_fmadd_pd(term2, x2_vec, C0_v); // ((...))*x2 + C0

    return _mm256_mul_pd(term1, x_vec); // x * poly_val
}

double calculate_s_avx() ;

// ======================================================================
// 测例 4: 使用 AVX 指令 + 4次循环展开
// ======================================================================

double calculate_s_avx_unrolled_x4();
