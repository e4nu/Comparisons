// Leave this at the top to enable features detected at build time in headers in
// HepMC3
#include <TH1D.h>
#include "NuHepMC/HepMC3Features.hxx"
#include "NuHepMC/EventUtils.hxx"
#include "NuHepMC/ReaderUtils.hxx"
#include "NuHepMC/WriterUtils.hxx"
#include "NuHepMC/make_writer.hxx"
#include "HepMC3/ReaderFactory.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "RadiativeCorrUtils.h"
#include "Utils.h"

///////////////////////////////////////////////////////////////////////////////
// comparison_1p1pi.cxx                                                      //
// --input-hepmc3-file : input MC file in hepmc3 format                      //
// --output-file : output MC file name, without format type. Def:myradevents //
///////////////////////////////////////////////////////////////////////////////

using namespace e4nu;
using namespace utils;

int main(int argc, char* argv[]) {

  std::string input_hepmc3_file = "";
  std::string output_name = "myradevents";
  int nevents = -1 ; // all
  double Delta_Em = 0.01;
  // process options
  if( argc > 1 ) { // configure rest of analysis
    if( utils::ExistArg("input-hepmc3-file",argc,argv)) {
      input_hepmc3_file = utils::GetArg("input-hepmc3-file",argc,argv);
    } else { std::cout << " --input-hepmc3-file is not defined "; return 0 ;}
    if( utils::ExistArg("output-file",argc,argv)) {
      output_name = utils::GetArg("output-file",argc,argv);
    }
  }

  auto rdr = HepMC3::deduce_reader(input_hepmc3_file);
  if (!rdr) {
    std::cout << "Failed to instantiate HepMC3::Reader from " << input_hepmc3_file << std::endl;
    return 1;
  }

  HepMC3::GenEvent evt;
  rdr->read_event(evt);
  if (rdr->failed()) {
    std::cout << "Failed to read first event from " << input_hepmc3_file << "."
              << std::endl;
    return 1;
  }

  auto in_gen_run_info = evt.run_info();
  auto vtx_statuses = NuHepMC::GR5::ReadVertexStatusIdDefinitions(in_gen_run_info);
  auto part_statuses = NuHepMC::GR6::ReadParticleStatusIdDefinitions(in_gen_run_info);
  auto out_gen_run_info = std::make_shared<HepMC3::GenRunInfo>(*in_gen_run_info);

  // re-open the file so that you start at the beginning
  rdr = HepMC3::deduce_reader(input_hepmc3_file);
  if (!rdr) {
    std::cout << "Failed to instantiate HepMC3::Reader from " << input_hepmc3_file << std::endl;
    return 1;
  }

  size_t nprocessed = 0;
  while (true) { // loop while there are events
    rdr->read_event(evt);
    if (rdr->failed()) {
      std::cout << "Reached the end of the file after " << nprocessed
                << " events." << std::endl;
      break;
    }
    evt.set_run_info(out_gen_run_info);
    evt.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);

    auto beampt = NuHepMC::Event::GetBeamParticle(evt);
    auto tgtpt = NuHepMC::Event::GetTargetParticle(evt);

    if (!beampt || !tgtpt) {  // this event didn't have a beam particle or target, its an odd one
      wrtr->write_event(evt); // write out events that we don't modify
      continue;
    }
    // Here anlaysis
    ++nprocessed;
    if( nevents > 0 && nprocessed > nevents ) break;
  }

  wrtr->close();
}
