from sklearn.metrics.pairwise import haversine_distances
from math import radians

# Example usage taken from docs:
# https://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise.haversine_distances.html
new_york = [40.7128, -74.0060]
boston = [42.3601, -71.0589]
new_york_in_radians = [radians(_) for _ in new_york]
boston_in_radians = [radians(_) for _ in boston]

# Compute the distance
for i in range(100_000):
    result = haversine_distances([new_york_in_radians, boston_in_radians])

# Multiply by Earth radius to get meters
distance = float(result[0][1]) * 6371000
print(f"Distance: {distance} meters")
