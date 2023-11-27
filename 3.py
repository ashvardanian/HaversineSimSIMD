import math
from numba import jit

EARTH_RADIUS = 6371000  # Earth radius in meters


@jit(nopython=True)
def to_radians(degrees):
    return math.radians(degrees)


@jit(nopython=True)
def haversine_distance(lat1, lon1, lat2, lon2):
    # Convert from degrees to radians
    lat1_rad = to_radians(lat1)
    lon1_rad = to_radians(lon1)
    lat2_rad = to_radians(lat2)
    lon2_rad = to_radians(lon2)

    # Difference in coordinates
    delta_lat = lat2_rad - lat1_rad
    delta_lon = lon2_rad - lon1_rad

    # Modified Haversine formula
    delta_lat_factor = math.sin(delta_lat / 2) ** 2
    mean_lat_factor = math.sin((lat1_rad + lat2_rad) / 2) ** 2
    delta_lon_factor = math.sin(delta_lon / 2) ** 2

    a = delta_lat_factor + (1 - delta_lat_factor - mean_lat_factor) * delta_lon_factor
    c = 2 * math.asin(math.sqrt(a))
    return EARTH_RADIUS * c


# Example usage
lat1 = 40.7128  # New York City latitude
lon1 = -74.0060  # New York City longitude
lat2 = 42.3601  # Boston latitude
lon2 = -71.0589  # Boston longitude

# Compute the distance
for i in range(10_000_000):
    distance = haversine_distance(lat1, lon1, lat2, lon2)
print(f"Distance: {distance} meters")
