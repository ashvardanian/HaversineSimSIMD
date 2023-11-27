import math

EARTH_RADIUS = 6371000  # Earth radius in meters


# Convert degrees to radians
def to_radians(degrees):
    return math.radians(degrees)


# Haversine formula to calculate distance
def haversine_distance(lat1, lon1, lat2, lon2):
    # Convert from degrees to radians
    lat1_rad = to_radians(lat1)
    lon1_rad = to_radians(lon1)
    lat2_rad = to_radians(lat2)
    lon2_rad = to_radians(lon2)

    # Difference in coordinates
    delta_lat = lat2_rad - lat1_rad
    delta_lon = lon2_rad - lon1_rad

    # Haversine formula
    a = (
        math.sin(delta_lat / 2) ** 2
        + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(delta_lon / 2) ** 2
    )
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # Distance
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
