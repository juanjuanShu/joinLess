#include "Common.h"

bool Common::hasRelation(const LocationType& location1, const LocationType& location2, double star_radius)
{
	return pow(location1.first - location2.first, 2) + pow(location1.second - location2.second, 2) <= star_radius * star_radius;
}
