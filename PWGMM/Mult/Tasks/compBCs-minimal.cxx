#include <Framework/AnalysisDataModel.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

// define columns in a sub-namespace of o2::aod
// and tables in namespace o2::aod
namespace o2::aod
{
namespace ambii
{
DECLARE_SOA_INDEX_COLUMN(MFTTrack, track);
DECLARE_SOA_INDEX_COLUMN(AmbiguousMFTTrack, ambMFTtrack);
//DECLARE_SOA_ARRAY_INDEX_COLUMN(AmbiguousMFTTrack, ambMFTtracks);
} // namespace ambii
DECLARE_SOA_INDEX_TABLE_USER(MA2T, MFTTracks, "MA2T", ambii::MFTTrackId, ambii::AmbiguousMFTTrackId);
//DECLARE_SOA_TABLE(MA2T, "AOD", "MA2T", ambii::AmbiguousMFTTrackIds);
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;

struct ProduceMFTAmbii {
  // build the index table MA2T
  Builds<aod::MA2T> idx;
  void init(InitContext const&){};
  // Produces<aod::MA2T> tambii;
  //
  // std::vector<int> iambii;
  // void process(aod::MFTTrack const& track, soa::SmallGroups<aod::AmbiguousMFTTracks> const& atracks)
  // {
  //   iambii.clear();
  //   for (auto& atrack : atracks) {
  //     iambii.push_back(atrack.globalIndex());
  //
  //   }
  //   tambii(iambii);
  // }
};

struct PrintMFTAmbii {
  void process(soa::Join<aod::MFTTracks,aod::MA2T> const& tracks, aod::AmbiguousMFTTracks const&)
  {
    for (auto& track : tracks) {
      if (track.has_ambMFTtrack()) {
        auto atrack = track.ambMFTtrack();
        printf("ambiguous mft track ID %lld, mft track ID %d\n", track.globalIndex(), atrack.mfttrackId());


        //LOGP(info, "MFT Track {} has {} entries", track.globalIndex(), track.ambMFTtracks().size());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ProduceMFTAmbii>(cfgc),
    adaptAnalysisTask<PrintMFTAmbii>(cfgc),
  };
}
