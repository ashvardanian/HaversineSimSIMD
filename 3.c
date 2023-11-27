#include <math.h>
#include <stdio.h>

float coords[2][2] = {
    {40.7128, -74.0060}, // New York City
    {42.3601, -71.0589}, // Boston
};

#define PI 3.14159265358979323846
#define EARTH_RADIUS 6371000 // in meters

// Convert degrees to radians
float to_radians(float degrees) { return degrees * (PI / 180.0); }

// Haversine formula to calculate distance
inline float haversine(float lat1, float lon1, float lat2, float lon2) {

  // Convert from degrees to radians
  float lat1_rad = to_radians(lat1);
  float lon1_rad = to_radians(lon1);
  float lat2_rad = to_radians(lat2);
  float lon2_rad = to_radians(lon2);

  // Difference in coordinates
  float delta_lat = lat2_rad - lat1_rad;
  float delta_lon = lon2_rad - lon1_rad;

  // Haversine formula
  float a = sinf(delta_lat / 2) * sinf(delta_lat / 2) +
            cosf(lat1_rad) * cosf(lat2_rad) * sinf(delta_lon / 2) *
                sinf(delta_lon / 2);

  float c = 2 * asinf(sqrtf(a));
  return EARTH_RADIUS * c;
}

inline float haversine_inequality(float lat1, float lon1, float lat2,
                                  float lon2) {

  // Convert from degrees to radians
  float lat1_rad = to_radians(lat1);
  float lon1_rad = to_radians(lon1);
  float lat2_rad = to_radians(lat2);
  float lon2_rad = to_radians(lon2);

  // Difference in coordinates
  float delta_lat = lat2_rad - lat1_rad;
  float delta_lon = lon2_rad - lon1_rad;

  // Haversine formula
  float a = sinf(delta_lat / 2) * sinf(delta_lat / 2) +
            cosf(lat1_rad) * cosf(lat2_rad) * sinf(delta_lon / 2) *
                sinf(delta_lon / 2);

  return a;
}

int main() {
  float distance;
  for (int it = 0; it < 100000000; it++) {
    float lat1 = coords[it % 2][0];
    float lon1 = coords[it % 2][1];
    float lat2 = coords[(it + 1) % 2][0];
    float lon2 = coords[(it + 1) % 2][1];
    distance = haversine(lat1, lon1, lat2, lon2);
    asm volatile("" : "+r"(distance)::);
  }

  printf("Distance: %f meters\n", distance);
  return 0;
}
