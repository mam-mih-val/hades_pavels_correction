//
// Created by mikhail on 11/30/20.
//

#ifndef HADES_RAPIDITY_SRC_RAPIDITY_H_
#define HADES_RAPIDITY_SRC_RAPIDITY_H_

#include <memory>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TProfile2D.h>
#include <TProfile3D.h>

#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/Matching.hpp>

#include <at_task/Task.h>
#include <DataContainer.hpp>
#include <QVector.hpp>

class Yield : public UserFillTask {

public:
  void UserInit(std::map<std::string, void *> &Map) override;
  void UserExec() override;
  void UserFinish() override;
  boost::program_options::options_description GetBoostOptions() override;
  void PreInit() override;
  void PostFinish() override {
    UserTask::PostFinish();
  }
private:
  double AngleDifference( double phi, double psi ){
    auto delta_phi = phi - psi;
    if( delta_phi < -M_PI )
      delta_phi+=2*M_PI;
    if( delta_phi > M_PI )
      delta_phi-=2*M_PI;
    return delta_phi;
  }
  static double ConeSideSquare(double theta1, double theta2, double a=100.0, double alpha=55.0*M_PI/180.0){
    double S=0;
    S = M_PI* a*a * cos(alpha)*cos(alpha) / sin(alpha) *
        (sin(theta1) / sin( alpha+theta1 )*sin(theta1) / sin( alpha+theta1 ) -
         sin(theta2) / sin( alpha+theta2 )*sin(theta2) / sin( alpha+theta2 ) );
    return S;
  }
  void LoopRecTracks();
  void LoopTruParticles();

  ATI2::Branch* event_header_{nullptr};
  ATI2::Branch* sim_header_{nullptr};
  ATI2::Branch* tracks_{nullptr};
  ATI2::Branch* sim_particles_{nullptr};

  AnalysisTree::Matching* rec_sim_matching_{nullptr};

  TH1F* h1_centrality_;

  // Average event shape of all particles
  TH3F* h3_tru_delta_phi_theta_centrality_all_{nullptr};
  TH3F* h3_rec_delta_phi_theta_centrality_all_{nullptr};
  // Average event shape of particles of interest
  TH3F* h3_tru_delta_phi_theta_centrality_pid_{nullptr};
  TH3F* h3_rec_delta_phi_theta_centrality_pid_{nullptr};
  // Additional QA distributions
  TH3F* h3_tru_pT_theta_centrality_pid_{nullptr};
  TH3F* h3_rec_pT_theta_centrality_pid_{nullptr};
  // Flow of particles
  TProfile2D* p2_tru_v1_pid_{nullptr};
  TProfile2D* p2_rec_v1_pid_{nullptr};
  TProfile2D* p2_tru_v1_all_{nullptr};
  TProfile2D* p2_rec_v1_all_{nullptr};
  // QA histograms
  TH1F* h1_psi_ep_;
  TH1F* h1_qx_;
  TH1F* h1_qy_;
  // MC QA histogram
  TH1F* h1_psi_rp_difference_;

  std::string str_qvector_file_name_;
  TFile* file_qvector_;
  TTree* tree_qvector_;
  std::string str_qvector_name_;
  Qn::DataContainer<Qn::QVector>* dc_qvector_;

  int reference_pdg_code_;

  bool is_mc_{false};

  int qvector_event_{0};
  int current_event_{0};

TASK_DEF(Yield, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
