#include <benchmark/benchmark.h>
#include <math.h>
#include <stdio.h>

#define PI 3.14159265358979323846
#define EARTH_RADIUS 6371000 // in meters

float coords[2][2]{
    {40.7128, -74.0060}, // New York City
    {42.3601, -71.0589}, // Boston
};

// Convert degrees to radians
float to_radians(float degrees) { return degrees * (PI / 180.0); }

inline float haversine_libc(float lat1, float lon1, float lat2, float lon2)
{
  // Convert from degrees to radians
  float lat1_rad = to_radians(lat1);
  float lon1_rad = to_radians(lon1);
  float lat2_rad = to_radians(lat2);
  float lon2_rad = to_radians(lon2);

  // Difference in coordinates
  float delta_lat = lat2_rad - lat1_rad;
  float delta_lon = lon2_rad - lon1_rad;

  // Haversine formula
  float a =
      sinf(delta_lat / 2) * sinf(delta_lat / 2) +
      cosf(lat1_rad) * cosf(lat2_rad) * sinf(delta_lon / 2) * sinf(delta_lon / 2);

  // This function is monotonically rising
  // https://www.wolframalpha.com/input?i=atan2%28sqrt%28x%29%2Csqrt%281-x%29%29
  float c = 2 * atan2f(sqrtf(a), sqrtf(1 - a));

  return EARTH_RADIUS * c;
}

inline float haversine_inequality(float lat1, float lon1, float lat2,
                                  float lon2)
{
  // Convert from degrees to radians
  float lat1_rad = to_radians(lat1);
  float lon1_rad = to_radians(lon1);
  float lat2_rad = to_radians(lat2);
  float lon2_rad = to_radians(lon2);

  // Difference in coordinates
  float delta_lat = lat2_rad - lat1_rad;
  float delta_lon = lon2_rad - lon1_rad;

  // Haversine formula
  float a =
      sinf(delta_lat / 2) * sinf(delta_lat / 2) +
      cosf(lat1Rad) * cosf(lat2_rad) * sinf(delta_lon / 2) * sinf(delta_lon / 2);

  return a;
}

static void libc(benchmark::State &state)
{
  size_t it = 0;
  for (auto _ : state)
  {
    float lat1 = coords[it % 2][0];
    float lon1 = coords[it % 2][1];
    float lat2 = coords[(it + 1) % 2][0];
    float lon2 = coords[(it + 1) % 2][1];
    benchmark::DoNotOptimize(haversine_libc(lat1, lon1, lat2, lon2));
    ++it;
  }
}
BENCHMARK(libc);

static void inequality(benchmark::State &state)
{
  size_t it = 0;
  for (auto _ : state)
  {
    float lat1 = coords[it % 2][0];
    float lon1 = coords[it % 2][1];
    float lat2 = coords[(it + 1) % 2][0];
    float lon2 = coords[(it + 1) % 2][1];
    benchmark::DoNotOptimize(haversine_inequality(lat1, lon1, lat2, lon2));
    ++it;
  }
}
BENCHMARK(inequality);

#ifdef __aarch64__

#include <arm_neon.h>

float haversine_simd(float lat1, float lon1, float lat2, float lon2)
{

  union
  {
    float32x4_t v;
    float f32s[4];
  } data;

  // Load and convert from degrees to radians.
  data.f32s[0] = lat1;
  data.f32s[1] = lon1;
  data.f32s[2] = lat2;
  data.f32s[3] = lon2;
  data.v = vmulq_n_f32(data.v, PI / 180.0);

  // Now we want to go from {lat1, lon1, lat2, lon2} to three values:
  //    ( ( lat2 - lat1 ) / 2) ^ 2
  //    ( ( lon2 - lon1 ) / 2) ^ 2
  //    ( ( lat2 + lat1 ) / 2) ^ 2
  // We are not gonna need the original values anymore, we can overwrite them.
  data.f32s[2] -= data.f32s[0];                // ( lat2 - lat1 )
  data.f32s[3] -= data.f32s[1];                // ( lon2 - lon1 )
  data.f32s[0] += data.f32s[2] + data.f32s[0]; // ( lat2 + lat1 )

  // Divide all of them by 2.
  data.v = vmulq_n_f32(data.v, 0.5);

  // Now map every value to its sinf using Horner's method.
  // For `sinf(x)` the algorithm with 4 iterations may look like this:
  //    double s = 1;
  //    s = 1 - s * ((x * x) / ((2 * 4 - 1) * (2 * 4 - 2)));
  //    s = 1 - s * ((x * x) / ((2 * 3 - 1) * (2 * 3 - 2)));
  //    s = 1 - s * ((x * x) / ((2 * 2 - 1) * (2 * 2 - 2)));
  //    s = s * x;
  //    x = s;
  // Precompute the constants:
  //    1 / ((2 * 4 - 1) * (2 * 4 - 2)) = 0.023809523809523808
  //    1 / ((2 * 3 - 1) * (2 * 3 - 2)) = 0.05
  //    1 / ((2 * 2 - 1) * (2 * 2 - 2)) = 0.16666666666666666

  float32x4_t x_squared = vmulq_f32(data.v, data.v);
  float32x4_t s = vdupq_n_f32(1);
  s = vsubq_f32(vdupq_n_f32(1),
                vmulq_n_f32(vmulq_f32(x_squared, s), 0.023809523809523808));
  s = vsubq_f32(vdupq_n_f32(1), vmulq_n_f32(vmulq_f32(x_squared, s), 0.05));
  s = vsubq_f32(vdupq_n_f32(1),
                vmulq_n_f32(vmulq_f32(x_squared, s), 0.16666666666666666));
  s = vmulq_f32(s, data.v);
  data.v = s;

  // Square the sines.
  data.v = vmulq_f32(data.v, data.v);

  // Final scalar reduction.
  float a = data.f32s[2] + (1 - data.f32s[2] - data.f32s[0]) * data.f32s[3];
  // float c = 2 * atan2f(sqrtf(a), sqrtf(1 - a));
  // return EARTH_RADIUS * c;
  return a;
}

#else

#include <xmmintrin.h> // SSE intrinsics

float haversine_simd(float lat1, float lon1, float lat2, float lon2)
{
  // Data union for alignment and ease of access
  union
  {
    __m128 v;
    float f32s[4];
  } data;

  // Load and convert from degrees to radians.
  data.f32s[0] = lat1;
  data.f32s[1] = lon1;
  data.f32s[2] = lat2;
  data.f32s[3] = lon2;
  data.v = _mm_mul_ps(data.v, _mm_set1_ps(PI / 180.0f));

  // Compute differences and sum
  data.f32s[2] -= data.f32s[0];                // ( lat2 - lat1 )
  data.f32s[3] -= data.f32s[1];                // ( lon2 - lon1 )
  data.f32s[0] += data.f32s[2] + data.f32s[0]; // ( lat2 + lat1 )

  // Divide all of them by 2.
  data.v = _mm_mul_ps(data.v, _mm_set1_ps(0.5f));

  // Apply Horner's method for sinf approximation
  __m128 x_squared = _mm_mul_ps(data.v, data.v);
  __m128 s = _mm_set1_ps(1.0f);
  s = _mm_sub_ps(
      _mm_set1_ps(1.0f),
      _mm_mul_ps(_mm_mul_ps(x_squared, s), _mm_set1_ps(0.023809523809523808f)));
  s = _mm_sub_ps(_mm_set1_ps(1.0f),
                 _mm_mul_ps(_mm_mul_ps(x_squared, s), _mm_set1_ps(0.05f)));
  s = _mm_sub_ps(
      _mm_set1_ps(1.0f),
      _mm_mul_ps(_mm_mul_ps(x_squared, s), _mm_set1_ps(0.16666666666666666f)));
  s = _mm_mul_ps(s, data.v);
  data.v = s;

  // Square the sines.
  data.v = _mm_mul_ps(data.v, data.v);

  // Final scalar reduction.
  float a = data.f32s[2] + (1 - data.f32s[2] - data.f32s[0]) * data.f32s[3];

  // Return the computed value.
  return a;
}

#endif

static void simd(benchmark::State &state)
{
  size_t it = 0;
  for (auto _ : state)
  {
    float lat1 = coords[it % 2][0];
    float lon1 = coords[it % 2][1];
    float lat2 = coords[(it + 1) % 2][0];
    float lon2 = coords[(it + 1) % 2][1];
    benchmark::DoNotOptimize(haversine_simd(lat1, lon1, lat2, lon2));
    ++it;
  }
}
BENCHMARK(simd);

BENCHMARK_MAIN();