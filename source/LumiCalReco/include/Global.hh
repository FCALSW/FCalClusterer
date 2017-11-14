#ifndef Global_hh
#define Global_hh 1

#include <map>
#include <memory>
#include <vector>

class LumiCalHit;
class LCCluster;
class ProjectionInfo;
class VirtualCluster;

using CalHit    = std::shared_ptr<LumiCalHit>;
using VecCalHit = std::vector<CalHit>;
using VDouble   = std::vector<double>;
using VInt      = std::vector<int>;

using MapIntCalHit    = std::map<int, CalHit>;
using MapIntLCCluster = std::map<int, LCCluster>;

using MapIntVCalHit = std::map<int, VecCalHit>;
using MapIntVDouble = std::map<int, VDouble>;
using MapIntVInt    = std::map<int, VInt>;

using MapIntMapIntLCCluster = std::map<int, MapIntLCCluster>;
using MapIntMapIntVCalHit   = std::map<int, MapIntVCalHit>;
using MapIntMapIntVDouble   = std::map<int, MapIntVDouble>;
using MapIntMapIntVInt      = std::map<int, MapIntVInt>;
using MapIntProjectionInfo  = std::map<int, ProjectionInfo>;
using MapIntVirtualCluster  = std::map<int, VirtualCluster>;
using MapIntDouble          = std::map<int, double>;
using MapIntInt             = std::map<int, int>;

using MapIntMapIntCalHit = std::map<int, MapIntCalHit>;

using VMapIntCalHit         = std::vector<MapIntCalHit>;
using VMapIntInt            = std::vector<MapIntInt>;
using VMapIntLCCluster      = std::vector<MapIntLCCluster>;
using VMapIntVInt           = std::vector<MapIntVInt>;
using VMapIntVirtualCluster = std::vector<MapIntVirtualCluster>;

using VVDouble = std::vector<VDouble>;

#define _VERY_SMALL_NUMBER 1e-9
#define _CLUSTER_MERGE_DEBUG 0
#define _GLOBAL_COUNTERS_UPDATE_DEBUG 0
#define _BHABHA_SELECTION_CUTS_DEBUG 1
#define _CLUSTER_STORE_DEBUG 0
#define _BP_DEBUG 0

#endif  // Global_hh
