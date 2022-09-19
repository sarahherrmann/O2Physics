#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace track
{
DECLARE_SOA_INDEX_COLUMN_FULL(BestCollision, bestCollision, int32_t, Collisions, "");
DECLARE_SOA_COLUMN(BestDCAXY, bestDCAXY, float);
DECLARE_SOA_COLUMN(BestDCAZ, bestDCAZ, float);
DECLARE_SOA_COLUMN(PtStatic, pts, float);
DECLARE_SOA_COLUMN(PStatic, ps, float);
DECLARE_SOA_COLUMN(EtaStatic, etas, float);
DECLARE_SOA_COLUMN(PhiStatic, phis, float);
} // namespace track
namespace fwdtrack
{
DECLARE_SOA_INDEX_COLUMN_FULL(BestCollision, bestCollision, int32_t, Collisions, "");
DECLARE_SOA_COLUMN(BestDCAXY, bestDCAXY, float);
DECLARE_SOA_COLUMN(PtStatic, pts, float);
DECLARE_SOA_COLUMN(PStatic, ps, float);
DECLARE_SOA_COLUMN(EtaStatic, etas, float);
DECLARE_SOA_COLUMN(PhiStatic, phis, float);
} // namespace fwdtrack
DECLARE_SOA_TABLE(BestCollisions, "AOD", "BESTCOLL",
                  aod::track::BestCollisionId, aod::track::BestDCAXY,
                  aod::track::BestDCAZ, track::X, track::Alpha, track::Y,
                  track::Z, track::Snp, track::Tgl, track::Signed1Pt,
                  track::PtStatic, track::PStatic, track::EtaStatic,
                  track::PhiStatic);

DECLARE_SOA_TABLE(BestCollisionsFwd, "AOD", "BESTCOLLFWD",
                  aod::fwdtrack::BestCollisionId, aod::fwdtrack::BestDCAXY,
                  fwdtrack::X, fwdtrack::Y,
                  fwdtrack::Z, fwdtrack::Tgl, fwdtrack::Signed1Pt,
                  fwdtrack::PtStatic, fwdtrack::PStatic, fwdtrack::EtaStatic,
                  fwdtrack::PhiStatic); // Snp does not exist
namespace indices
{
DECLARE_SOA_ARRAY_INDEX_COLUMN(Collision, collisions);
}

DECLARE_SOA_TABLE(MatchedMulti, "AOD", "MAMU",
                  indices::BCId, indices::CollisionIds);
} // namespace o2::aod
