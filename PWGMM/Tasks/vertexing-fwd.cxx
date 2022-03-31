#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "MathUtils/Utils.h"
#include "CommonConstants/LHCConstants.h"

using namespace o2;
using namespace o2::framework;

struct vertexingfwd {
  

  // std::vector<double> propagateToZlinear(double zEnd, aod::AmbiguousTrack const& track)
  // {
  //   // Track parameters and their covariances linearly extrapolated to the plane at "zEnd".
  //
  //   //Calculate the jacobian related to the track parameters extrapolated to "zEnd"
  //   auto dZ = (zEnd - track.z());
  //   auto phi0 = track.phi();
  //   auto tanl0 = track.tgl();
  //   auto invtanl0 = 1.0 / tanl0;
  //   auto [sinphi0, cosphi0] = o2::math_utils::sincosd(phi0);
  //   auto n = dZ * invtanl0;
  //
  //   std::vector<double> T ={track.x(),track.y(),zEnd};
  //   // Extrapolate track parameters to "zEnd"
  //   T[0] += n * cosphi0;
  //   T[1] += n * sinphi0;
  //
  //   //Update to do : Add the jacobian and turn it to a method if possible
  //   return(T);
  //
  // }

  // void process(aod::AmbiguousTracks const& tracks, aod::BCs const&, aod::Collisions const& collisions)//AmbiguousMFTTracks and fwd doesn't work
  // {
  //   std::vector<double> T;
  //   for (auto& track : tracks) {
  //     nCounter++;
  //
  //     LOGF(info, "------------------------------------ We look at track %d which has %d possible BCs", track.globalIndex(), track.bc().size());
  //     for (auto& bc : track.bc())
  //     {
  //       //LOGF(info, "  BC %d with global BC %lld", bc.globalIndex(), bc.globalBC());
  //       if (nCounter < nMax)
  //       {
  //         for (auto& collision : collisions)
  //         {
  //           //printf("collision BC ID = %lld, bc.globalIndex() %lld\n", collision.bcId(), bc.globalIndex());
  //           if (collision.bcId() == bc.globalIndex())
  //           {
  //             printf("collision BC ID = %lld, bc.globalIndex() %lld\n", collision.bcId(), bc.globalIndex());
  //             printf("collision pos Z %f\n", collision.posZ());
  //              T = propagateToZlinear(collision.posZ(), track);
  //             printf("T[0] %f, T[1] %f\n", T[0], T[1]);
  //           }
  //         }
  //       }
  //     }
  //
  //   }
  // }
  void process(aod::AmbiguousTracks const& ambitracks, aod::Tracks const& tracks)//AmbiguousMFTTracks and fwd doesn't work
  {
    for (auto& ambitrack : ambitracks)
    {
      auto extAmbiTrack =tracks.iteratorAt(ambitrack.globalIndex());
      printf("phi = %f\n", extAmbiTrack.phi());
    }


  }


};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<vertexingfwd>(cfgc)};
}
