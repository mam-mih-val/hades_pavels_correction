//
// Created by mikhail on 11/30/20.
//

#include "yield.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"

#include <AnalysisTree/DataHeader.hpp>

TASK_IMPL(Yield)

boost::program_options::options_description Yield::GetBoostOptions() {
  using namespace boost::program_options;
  options_description desc(GetName() + " options");
  desc.add_options()
    ("q-vector-file", value(&str_qvector_file_name_), "Path to the file containing Q-vectors")
    ("q-vector-name", value(&str_qvector_name_), "Path to the file containing Q-vectors")
    ("pdg-code", value(&reference_pdg_code_)->default_value(2212), "PDG-code of particle");
  return desc;
}

void Yield::PreInit() {
//  SetInputBranchNames({"mdc_vtx_tracks", "event_header", "sim_tracks"});
}

void Yield::UserInit(std::map<std::string, void *> &Map) {
  // extracting event variables
  event_header_ = GetInBranch("event_header");
  tracks_ = GetInBranch("mdc_vtx_tracks");
  try {
    sim_header_ = GetInBranch("sim_header");
    sim_particles_ = GetInBranch("sim_tracks");
    rec_sim_matching_ = static_cast<AnalysisTree::Matching *>(
        Map.at("mdc_vtx_tracks2sim_tracks"));
    is_mc_ = true;
  } catch (std::exception &) {
    is_mc_ = false;
  }

  h1_centrality_ = new TH1F( "centrality", ";TOF+RPC hits centrality (%)", 20, 0.0, 100.0 );

  h3_rec_delta_phi_theta_centrality_all_ = new TH3F("h3_rec_delta_phi_theta_centrality_all",
                                                    ";#Delta#phi (rad);#theta;centrality (%)",
                                                    350, -3.5, 3.5,
                                                    140, 0.2, 1.6,
                                                    12, 0.0, 60.);
  h3_rec_delta_phi_theta_centrality_pid_ = new TH3F("h3_rec_delta_phi_theta_centrality_pid",
                                                    ";#Delta#phi (rad);#theta;centrality (%)",
                                                    350, -3.5, 3.5,
                                                    140, 0.2, 1.6,
                                                    12, 0.0, 60.);

  h3_rec_pT_theta_centrality_pid_ = new TH3F( "h3_rec_pT_theta_centrality_pid",
                                          ";p_{T} (GeV/c);#theta (rad);centrality (%)",
                                          200, 0.0, 2.0,
                                          170, 0.0, 1.7,
                                          12, 0.0, 60.);

  p2_rec_v1_pid_ = new TProfile2D( "p2_rec_v1_pid", ";theta;centrality", 140, 0.2, 1.6, 12, 0.0, 60 );

  p2_rec_v1_all_ = new TProfile2D( "p2_rec_v1_all", ";theta;centrality", 140, 0.2, 1.6, 12, 0.0, 60 );

  if( is_mc_ ) {
    h3_tru_delta_phi_theta_centrality_all_ =
        new TH3F("h3_tru_delta_phi_theta_centrality_all",
                 ";#Delta#phi (rad);#theta;centrality (%)", 350, -3.5, 3.5, 140,
                 0.2, 1.6, 12, 0.0, 60.);
    h3_tru_delta_phi_theta_centrality_pid_ =
        new TH3F("h3_tru_delta_phi_theta_centrality_pid",
                 ";#Delta#phi (rad);#theta;centrality (%)", 350, -3.5, 3.5, 140,
                 0.2, 1.6, 12, 0.0, 60.);
    h3_tru_pT_theta_centrality_pid_ =
        new TH3F("h3_tru_pT_theta_centrality_pid",
                 ";p_{T} (GeV/c);#theta (rad);centrality (%)", 200, 0.0, 2.0,
                 170, 0.0, 1.7, 12, 0.0, 60.);
    p2_tru_v1_pid_ = new TProfile2D("p2_tru_v1_pid", ";theta;centrality", 140,
                                    0.2, 1.6, 12, 0.0, 60);
    p2_tru_v1_all_ = new TProfile2D("p2_tru_v1_all", ";theta;centrality", 140,
                                    0.2, 1.6, 12, 0.0, 60);
  }

  file_qvector_ = TFile::Open(str_qvector_file_name_.c_str());
  if( !file_qvector_ )
    throw std::runtime_error( "File "+str_qvector_file_name_+" does not exist" );
  file_qvector_->GetObject("tree", tree_qvector_);
  if( !tree_qvector_ )
    throw std::runtime_error( "There is no tree in file "+str_qvector_file_name_ );
  tree_qvector_->SetBranchAddress(str_qvector_name_.c_str(), &dc_qvector_);
  if( !dc_qvector_ )
    throw std::runtime_error( "There is no branch "+str_qvector_name_+" in tree" );
  out_file_->cd();
  std::cout << "Initialized" << std::endl;
}

void Yield::UserExec() {
  if( current_event_ == 0 ) {
    current_event_++;
    return;
  }
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  h1_centrality_->Fill( centrality );
  tree_qvector_->GetEntry(qvector_event_);
  this->LoopRecTracks();
  if( is_mc_ )
    this->LoopTruParticles();
  qvector_event_++;
  current_event_++;
}

void Yield::LoopRecTracks() {
  using AnalysisTree::Particle;
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();

  auto rec_chi2_var = GetVar("mdc_vtx_tracks/chi2");
  auto rec_dca_xy_var = GetVar("mdc_vtx_tracks/dca_xy");
  auto rec_dca_z_var = GetVar("mdc_vtx_tracks/dca_z");

  auto qvec = dc_qvector_->At(0);
  auto psi_rp = qvec.psi(1);

  for (auto track : tracks_->Loop()) {
    auto pid = track.DataT<Particle>()->GetPid();
    auto mass = track.DataT<Particle>()->GetMass();
    if( pid != 0 ) {
      if( TDatabasePDG::Instance()->GetParticle(pid) )
        mass = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
    }
    auto mom4 = track.DataT<Particle>()->Get4MomentumByMass(mass);
    auto chi2 = track[rec_chi2_var].GetVal();
    auto dca_xy = track[rec_dca_xy_var].GetVal();
    auto dca_z = track[rec_dca_z_var].GetVal();
    auto delta_phi = AngleDifference(mom4.Phi(), psi_rp);
    if( mom4.Pt() > 0.4 ) {
      h3_rec_delta_phi_theta_centrality_all_->Fill(delta_phi, mom4.Theta(),
                                                   centrality);
    }
    p2_rec_v1_all_->Fill( mom4.Theta(), centrality, cos(delta_phi) );
    if( chi2 > 100.0 )
      continue;
    if ( -10 > dca_xy || dca_xy > 10 )
      continue;
    if ( -10 > dca_z || dca_z > 10 )
      continue;
    if( pid != reference_pdg_code_ )
      continue;
    h3_rec_pT_theta_centrality_pid_->Fill(mom4.Pt(), mom4.Theta(), centrality);
    if( pid == 2212 && mom4.Pt() < 0.4 )
      continue;
    if( abs(pid) == 211 && mom4.Pt() < 0.2 )
      continue;
    h3_rec_delta_phi_theta_centrality_pid_->Fill(delta_phi, mom4.Theta(), centrality);
    p2_rec_v1_pid_->Fill( mom4.Theta(), centrality, cos(delta_phi) );
  }
}

void Yield::LoopTruParticles() {
  using AnalysisTree::Particle;
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();

  auto var_is_primary = GetVar("sim_tracks/is_primary");
  auto qvec = dc_qvector_->At(0);
  auto psi_rp = qvec.psi(1);

  int idx1 =-1;
  for( auto particle : sim_particles_->Loop() ){
    idx1++;
    auto mass = particle.DataT<Particle>()->GetMass();
    auto pid = particle.DataT<Particle>()->GetPid();
    auto is_prim = particle[var_is_primary].GetBool();
    auto mom4 = particle.DataT<Particle>()->Get4MomentumByMass(mass);
    auto phi = mom4.Phi() + M_PI;
    double charge=0.0;
    if( TDatabasePDG::Instance()->GetParticle( pid ) ){
      charge= TDatabasePDG::Instance()->GetParticle( pid )->Charge() / 3.0;
    }
    auto delta_phi = AngleDifference(mom4.Phi(), psi_rp);
    if( fabs(charge) < 0.01 )
      continue;
    if( mom4.Pt() > 0.4 ) {
      h3_tru_delta_phi_theta_centrality_all_->Fill(delta_phi, mom4.Theta(),
                                                   centrality);
    }
    p2_tru_v1_all_->Fill( mom4.Theta(), centrality, cos(delta_phi) );
    if( !is_prim )
      continue;
    if( pid!=reference_pdg_code_ )
      continue;
    h3_tru_pT_theta_centrality_pid_->Fill(mom4.Pt(), mom4.Theta(), centrality);
    if( pid == 2212 && mom4.Pt() < 0.4 )
      continue;
    if( abs(pid) == 211 && mom4.Pt() < 0.2 )
      continue;
    h3_tru_delta_phi_theta_centrality_pid_->Fill(delta_phi, mom4.Theta(), centrality);
    p2_tru_v1_pid_->Fill( mom4.Theta(), centrality, cos(delta_phi) );
  }
}


void Yield::UserFinish() {
  out_file_->cd();
  h1_centrality_->Write();

//  out_file_->mkdir("efficiency_projections");
//  out_file_->cd("efficiency_projections");

  out_file_->cd();
  h3_rec_delta_phi_theta_centrality_all_->Write();
  h3_rec_delta_phi_theta_centrality_pid_->Write();
  h3_rec_pT_theta_centrality_pid_->Write();
  p2_rec_v1_pid_->Write();
  p2_rec_v1_all_->Write();
  if( is_mc_ ) {
    h3_tru_delta_phi_theta_centrality_all_->Write();
    h3_tru_delta_phi_theta_centrality_pid_->Write();
    h3_tru_pT_theta_centrality_pid_->Write();
    p2_tru_v1_pid_->Write();
    p2_tru_v1_all_->Write();
  }
  std::cout << "Finished" << std::endl;
}
